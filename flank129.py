#!/usr/bin/env python3
"""flank129: read the 129-strain DNA around a targeted gene from a knockout study's own RNA-seq reads.

Most knockout alleles were made in 129 embryonic stem cells and backcrossed to C57BL/6, so knockouts carry 129 DNA around the targeted
gene that their controls may lack. flank129 genotypes every sample at the sites where the 129 strains 129S1/SvImJ, 129P2/OlaHsd and
129S5/SvEvBrd are all homozygous for the same non-reference allele (Mouse Genomes Project REL-1505, GRCm38), calls 1 Mb bins as
129-carrying, and reports the segments in which every knockout and no control sample carries 129 alleles (and the reverse).

Commands
  sites  VCF_129S1 VCF_129P2 VCF_129S5 OUT.tsv.gz     shared homozygous 129 sites from the three REL-1505 strain VCFs (needs bcftools)
  count  BAM SITES.tsv.gz OUT.ad.tsv.gz               reads carrying the C57BL/6 (reference) and the 129 base at every site (needs samtools;
                                                      BAM aligned to GRCm38, e.g. hisat2 -x grcm38/genome ... | samtools sort)
  call   SAMPLES.tsv GENES.gtf.gz TARGET OUTDIR       bins and segments for one study; SAMPLES.tsv has columns sample, group (knockout or
                                                      control) and file (one row per .ad.tsv.gz; several rows of one sample are summed)
         [--de DE.tsv]                                optional: differential expression of the study (columns gene, fdr): genes with
                                                      FDR < 0.05 inside the knockout-only segments near the target against all others

Rules (as in the paper): sites with at least 4 reads; a 1 Mb bin with at least 3 such sites is 129-carrying when its 129 allele fraction is
at least 0.2; a bin enters the comparison of groups only when every sample of both groups has a call there; adjacent separating bins are
merged; segments within 20 Mb of the target (or containing it) are reported as near the target.
"""
import argparse, gzip, multiprocessing, os, re, subprocess, sys
from concurrent.futures import ProcessPoolExecutor
import numpy as np, pandas as pd

MIN_READS, MIN_SITES, CALL, BIN, NEAR_MB = 4, 3, 0.2, 1_000_000, 20


def cmd_sites(a):
    sets = []
    for v in (a.vcf_129s1, a.vcf_129p2, a.vcf_129s5):
        q = subprocess.run(f"bcftools view -i 'FILTER=\"PASS\" && GT=\"1/1\"' -Ou {v} | bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n'",
                           shell=True, capture_output=True, text=True, check=True).stdout
        d = pd.DataFrame([l.split('\t') for l in q.splitlines()], columns=['chrom', 'pos', 'ref', 'alt']); d['pos'] = d.pos.astype(int); sets.append(d)
    s = sets[0].merge(sets[1], on=['chrom', 'pos', 'ref', 'alt']).merge(sets[2], on=['chrom', 'pos', 'ref', 'alt'])
    s = s[s.alt.str.len() == 1].sort_values(['chrom', 'pos']); s.to_csv(a.out, sep='\t', index=False)
    print(f'{len(s):,} shared homozygous 129 sites -> {a.out}')


CARET, INDEL = re.compile(r'\^.', re.S), re.compile(r'[+-](\d+)')


def _pile(args):
    """samtools mpileup over one chromosome at the sites (MAPQ >= 20, base quality >= 20, no BAQ); reads with the reference and the 129 base."""
    bam, chrom, bed, refalt, samtools = args; out = []
    pile = subprocess.run(f'{samtools} mpileup -r {chrom} -l {bed} -q 20 -Q 20 -d 100000 -B {bam} 2>/dev/null', shell=True, capture_output=True, text=True).stdout
    for line in pile.splitlines():
        c, pos, _, dp, bases = line.split('\t')[:5]
        ra = refalt.get((c, int(pos)))
        if ra is None or dp == '0': continue
        if '^' in bases: bases = CARET.sub('', bases)
        if '+' in bases or '-' in bases:
            parts, k = [], 0
            for m in INDEL.finditer(bases):
                if m.start() < k: continue
                parts.append(bases[k:m.start()]); k = m.end() + int(m.group(1))
            parts.append(bases[k:]); bases = ''.join(parts)
        b = bases.upper(); t = sum(b.count(x) for x in 'ACGTN'); nr, na = b.count(ra[0]), b.count(ra[1])
        out.append(f'{c}\t{pos}\t{ra[0]}\t{ra[1]}\t{nr}\t{na}\t{t - nr - na}\n')
    return ''.join(out)


def cmd_count(a):
    sites = pd.read_csv(a.sites, sep='\t', dtype={'chrom': str})
    bed = a.out + '.sites.bed'
    pd.DataFrame({'c': sites.chrom, 's': sites.pos - 1, 'e': sites.pos}).to_csv(bed, sep='\t', header=False, index=False)
    refalt = {(c, int(p)): (r, al) for c, p, r, al in zip(sites.chrom, sites.pos, sites.ref, sites.alt)}
    hdr = subprocess.run(f'{a.samtools} view -H {a.bam}', shell=True, capture_output=True, text=True, check=True).stdout
    chroms = [l.split('\t')[1][3:] for l in hdr.splitlines() if l.startswith('@SQ')]
    with ProcessPoolExecutor(a.threads, mp_context=multiprocessing.get_context('fork')) as ex:
        chunks = list(ex.map(_pile, [(a.bam, ch, bed, refalt, a.samtools) for ch in chroms]))
    with gzip.open(a.out, 'wt') as fo:
        fo.write('chrom\tpos\tref\talt\tn_ref\tn_alt\tn_other\n'); fo.write(''.join(chunks))
    os.remove(bed); print(f'{a.bam} -> {a.out}')


def bins_of(files):
    """1 Mb bins of one sample (all its files summed per site)."""
    x = pd.concat([pd.read_csv(f, sep='\t', usecols=['chrom', 'pos', 'n_ref', 'n_alt'], dtype={'chrom': str}) for f in files])
    x = x.groupby(['chrom', 'pos'], as_index=False)[['n_ref', 'n_alt']].sum()
    x = x[(x.n_ref + x.n_alt) >= MIN_READS].copy(); x['bin'] = x.pos // BIN; x['f'] = x.n_alt / (x.n_ref + x.n_alt)
    y = x.groupby(['chrom', 'bin']).agg(sites=('pos', 'size'), ref=('n_ref', 'sum'), alt=('n_alt', 'sum'), mean_site_frac=('f', 'mean')).reset_index()
    y = y[y.sites >= MIN_SITES].copy(); y['frac'] = y.alt / (y.ref + y.alt); y['is129'] = y.frac >= CALL
    return y


def cmd_call(a):
    S = pd.read_csv(a.samples, sep='\t'); S['group'] = S.group.str.lower()
    if not set(S.group) <= {'knockout', 'control'}: sys.exit('group must be knockout or control')
    g = pd.read_csv(a.gtf, sep='\t', header=None, comment='#', usecols=[0, 2, 3, 4, 8], names=['chrom', 'type', 'start', 'end', 'attr'], dtype={'chrom': str})
    g = g[g.type == 'gene']; g['symbol'] = g.attr.str.extract(r'gene_name "([^"]+)"'); g = g.drop(columns=['attr', 'type'])
    t = g[g.symbol == a.target]
    if t.empty: sys.exit(f'{a.target} not in {a.gtf}')
    t = t.iloc[0]; tmid = (t.start + t.end) / 2
    parts = []
    for (smp, grp), z in S.groupby(['sample', 'group'], sort=True):
        y = bins_of(list(z.file)); y.insert(0, 'group', grp); y.insert(0, 'sample', smp); parts.append(y)
    B = pd.concat(parts, ignore_index=True)
    nko, nct = S[S.group == 'knockout']['sample'].nunique(), S[S.group == 'control']['sample'].nunique()
    m = B.groupby(['chrom', 'bin', 'group']).agg(m=('is129', 'mean'), k=('sample', 'nunique')).unstack('group')
    m = m.reindex(columns=pd.MultiIndex.from_product([['m', 'k'], ['knockout', 'control']]))
    full = (m[('k', 'knockout')].fillna(0) == nko) & (m[('k', 'control')].fillna(0) == nct)
    s = pd.DataFrame({'knockout': m[('m', 'knockout')], 'control': m[('m', 'control')]})[full].reset_index()
    s['kind'] = np.where((s.knockout == 1) & (s.control == 0), 'knockout carries 129', np.where((s.knockout == 0) & (s.control == 1), 'control carries 129', ''))
    s = s[s.kind != ''].sort_values(['chrom', 'bin'])
    de = pd.read_csv(a.de, sep='\t') if a.de else None
    segs = []
    for (c, kind), z in s.groupby(['chrom', 'kind']):
        start = prev = None
        for bn in list(z.bin) + [None]:
            if bn is not None and prev is not None and bn == prev + 1: prev = bn; continue
            if start is not None:
                lo, hi = start * BIN, (prev + 1) * BIN
                dist = 0 if (c == t.chrom and lo <= tmid <= hi) else (min(abs(lo - tmid), abs(hi - tmid)) if c == t.chrom else np.inf)
                row = {'kind': kind, 'chrom': c, 'start_mb': lo / 1e6, 'end_mb': hi / 1e6, 'length_mb': (hi - lo) / 1e6, 'contains_target': dist == 0,
                       'distance_mb': dist / 1e6 if np.isfinite(dist) else None}
                if de is not None:
                    genes = g[(g.chrom == c) & (g.end >= lo) & (g.start <= hi)].symbol.dropna().unique(); dd = de[de.gene.isin(genes)]
                    row.update(genes_tested=len(dd), genes_fdr05=int((dd.fdr < 0.05).sum()))
                segs.append(row)
            start = prev = bn
    G = pd.DataFrame(segs, columns=['kind', 'chrom', 'start_mb', 'end_mb', 'length_mb', 'contains_target', 'distance_mb'] + (['genes_tested', 'genes_fdr05'] if de is not None else []))
    os.makedirs(a.outdir, exist_ok=True)
    B.to_csv(os.path.join(a.outdir, 'bins.tsv.gz'), sep='\t', index=False); G.to_csv(os.path.join(a.outdir, 'segments.tsv'), sep='\t', index=False)
    near = G[(G.chrom == t.chrom) & (G.distance_mb.fillna(np.inf) <= NEAR_MB)] if len(G) else G
    print(f'{a.target} (chromosome {t.chrom}, {tmid / 1e6:.2f} Mb); {nko} knockout and {nct} control samples; {B.groupby("sample").size().median():.0f} called bins per sample')
    for kind in ('knockout carries 129', 'control carries 129'):
        k = near[near.kind == kind] if len(near) else near
        print(f'  {kind} within {NEAR_MB} Mb of the target: {len(k)} segment(s), {k.length_mb.sum() if len(k) else 0:.0f} Mb'
              + (f'; contains the target' if len(k) and k.contains_target.any() else ''))
    if de is not None and len(near):
        k = near[near.kind == 'knockout carries 129']
        if len(k):
            ins = set()
            for r in k.itertuples(): ins |= set(g[(g.chrom == r.chrom) & (g.end >= r.start_mb * 1e6) & (g.start <= r.end_mb * 1e6)].symbol.dropna())
            dd = de.assign(inside=de.gene.isin(ins)); a_, b_ = dd[dd.inside], dd[~dd.inside]
            print(f'  differentially expressed (FDR < 0.05): {(a_.fdr < 0.05).sum()} of {len(a_)} genes inside the knockout-only segments, '
                  f'{(b_.fdr < 0.05).sum()} of {len(b_)} elsewhere')
    print(f'bins and segments -> {a.outdir}')


def main():
    ap = argparse.ArgumentParser(prog='flank129', description=__doc__.split('\n\n')[0]); sub = ap.add_subparsers(dest='cmd', required=True)
    p = sub.add_parser('sites'); p.add_argument('vcf_129s1'); p.add_argument('vcf_129p2'); p.add_argument('vcf_129s5'); p.add_argument('out'); p.set_defaults(f=cmd_sites)
    p = sub.add_parser('count'); p.add_argument('bam'); p.add_argument('sites'); p.add_argument('out'); p.add_argument('--threads', type=int, default=8)
    p.add_argument('--samtools', default='samtools'); p.set_defaults(f=cmd_count)
    p = sub.add_parser('call'); p.add_argument('samples'); p.add_argument('gtf'); p.add_argument('target'); p.add_argument('outdir'); p.add_argument('--de'); p.set_defaults(f=cmd_call)
    a = ap.parse_args(); a.f(a)


if __name__ == '__main__':
    main()
