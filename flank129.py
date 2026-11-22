#!/usr/bin/env python3
"""flank129: read the 129-strain DNA around a targeted gene from a knockout study's own RNA-seq reads.

Most knockout alleles were made in 129 embryonic stem cells and backcrossed to C57BL/6, so knockouts carry 129 DNA around the targeted
gene that their controls may lack. flank129 genotypes every sample at the sites where the 129 strains 129S1/SvImJ, 129P2/OlaHsd and
129S5/SvEvBrd are all homozygous for the same non-reference allele (Mouse Genomes Project REL-1505, GRCm38), calls 1 Mb bins as
129-carrying, and reports the segments in which the knockouts but not the controls carry 129 alleles (and the reverse). With a site
list made from another strain's VCF (`sites`), the same steps read DNA of that strain, e.g. DBA/2J or CBA/J around alleles edited in
hybrid zygotes.

Commands
  sites  VCF [VCF ...] OUT.tsv.gz                    sites where all given REL-1505 strain VCFs are homozygous for the same alternative
                                                      SNV (the 129 list: the 129S1, 129P2 and 129S5 VCFs; another strain: its VCF; needs bcftools)
  count  BAM SITES.tsv.gz OUT.ad.tsv.gz               reads carrying the C57BL/6 (reference) and the 129 base at every site (needs samtools;
                                                      BAM aligned to GRCm38, e.g. hisat2 -x grcm38/genome ... | samtools sort)
  call   SAMPLES.tsv GENES.gtf.gz TARGET OUTDIR       bins and segments for one study; SAMPLES.tsv has columns sample, group (knockout or
                                                      control) and file (one row per .ad.tsv.gz; several rows of one sample are summed)
         [--de DE.tsv]                                optional: differential expression of the study (columns gene, fdr): genes with
                                                      FDR < 0.05 inside the knockout-only segments near the target against all others
         [--strict]                                   the rule of version 1.0 (every sample called; only directly adjacent bins merged)

Rules (version 1.1, as in the paper): sites with at least 4 reads; a 1 Mb bin with at least 3 such sites is 129-carrying when its 129
allele fraction is at least 0.2. A bin enters the comparison of groups when each group has calls in at least max(2, half) of its samples;
it separates the groups when every called knockout and no called control is 129-carrying (or the reverse). Separating bins of the same kind
are joined across gaps of up to 5 bins that did not enter the comparison (a bin without enough reads in some samples does not break a
segment; a bin that entered with any other state does). Segments within 20 Mb of the target (or containing it) are near the target.
Version 1.0 required calls in every sample and joined only directly adjacent bins (--strict).
"""
import argparse, gzip, multiprocessing, os, re, subprocess, sys
from concurrent.futures import ProcessPoolExecutor
import numpy as np, pandas as pd

MIN_READS, MIN_SITES, CALL, BIN, NEAR_MB, MAX_GAP = 4, 3, 0.2, 1_000_000, 20, 5
__version__ = '1.1.0'


def cmd_sites(a):
    *vcfs, out = a.files
    if not vcfs: sys.exit('sites: give at least one VCF and the output file')
    sets = []
    for v in vcfs:
        q = subprocess.run(f"bcftools view -i 'FILTER=\"PASS\" && GT=\"1/1\"' -Ou {v} | bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n'",
                           shell=True, capture_output=True, text=True, check=True).stdout
        d = pd.DataFrame([l.split('\t') for l in q.splitlines()], columns=['chrom', 'pos', 'ref', 'alt']); d['pos'] = d.pos.astype(int)
        d['alt'] = d.alt.str.split(',').str[0]; sets.append(d)
    s = sets[0]
    for d in sets[1:]: s = s.merge(d, on=['chrom', 'pos', 'ref', 'alt'])
    s = s[(s.ref.str.len() == 1) & (s.alt.str.len() == 1)].sort_values(['chrom', 'pos']); s.to_csv(out, sep='\t', index=False)
    print(f'{len(s):,} shared homozygous sites from {len(vcfs)} VCF(s) -> {out}')


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


def group_segments(B, nko, nct, strict=False):
    """Segments separating the groups. B: one row per sample and called bin (columns sample, group, chrom, bin, is129).
    Returns (chrom, kind, first bin, last bin, bridged Mb) tuples, kind 'knockout carries 129' or 'control carries 129'."""
    m = B.groupby(['chrom', 'bin', 'group']).agg(m=('is129', 'mean'), k=('sample', 'nunique')).unstack('group')
    m = m.reindex(columns=pd.MultiIndex.from_product([['m', 'k'], ['knockout', 'control']]))
    kk, kc = m[('k', 'knockout')].fillna(0), m[('k', 'control')].fillna(0)
    enter = (kk == nko) & (kc == nct) if strict else (kk >= max(2, -(-nko // 2))) & (kc >= max(2, -(-nct // 2)))
    s = pd.DataFrame({'knockout': m[('m', 'knockout')], 'control': m[('m', 'control')]})[enter].reset_index()
    s['kind'] = np.where((s.knockout == 1) & (s.control == 0), 'knockout carries 129', np.where((s.knockout == 0) & (s.control == 1), 'control carries 129', ''))
    out, gap = [], (0 if strict else MAX_GAP)
    for c, z in s.sort_values(['chrom', 'bin']).groupby('chrom', sort=True):
        cur = None                                              # [kind, first, last, bridged]
        for bn, kind in zip(z.bin, z.kind):                     # every bin that entered the comparison, in order
            if cur is not None and kind == cur[0] and bn - cur[2] - 1 <= gap:
                cur[3] += bn - cur[2] - 1; cur[2] = bn; continue
            if cur is not None: out.append((c, cur[0], cur[1], cur[2], cur[3])); cur = None
            if kind: cur = [kind, bn, bn, 0]
        if cur is not None: out.append((c, cur[0], cur[1], cur[2], cur[3]))
    return out


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
    de = pd.read_csv(a.de, sep='\t') if a.de else None
    segs = []
    for c, kind, b0, b1, bridged in group_segments(B, nko, nct, strict=a.strict):
        lo, hi = b0 * BIN, (b1 + 1) * BIN
        dist = 0 if (c == t.chrom and lo <= tmid <= hi) else (min(abs(lo - tmid), abs(hi - tmid)) if c == t.chrom else np.inf)
        row = {'kind': kind, 'chrom': c, 'start_mb': lo / 1e6, 'end_mb': hi / 1e6, 'length_mb': (hi - lo) / 1e6, 'bridged_mb': bridged,
               'contains_target': dist == 0, 'distance_mb': dist / 1e6 if np.isfinite(dist) else None}
        if de is not None:
            genes = g[(g.chrom == c) & (g.end >= lo) & (g.start <= hi)].symbol.dropna().unique(); dd = de[de.gene.isin(genes)]
            row.update(genes_tested=len(dd), genes_fdr05=int((dd.fdr < 0.05).sum()))
        segs.append(row)
    G = pd.DataFrame(segs, columns=['kind', 'chrom', 'start_mb', 'end_mb', 'length_mb', 'bridged_mb', 'contains_target', 'distance_mb'] + (['genes_tested', 'genes_fdr05'] if de is not None else []))
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
    ap = argparse.ArgumentParser(prog='flank129', description=__doc__.split('\n\n')[0]); ap.add_argument('--version', action='version', version=__version__); sub = ap.add_subparsers(dest='cmd', required=True)
    p = sub.add_parser('sites'); p.add_argument('files', nargs='+', help='VCF [VCF ...] OUT.tsv.gz'); p.set_defaults(f=cmd_sites)
    p = sub.add_parser('count'); p.add_argument('bam'); p.add_argument('sites'); p.add_argument('out'); p.add_argument('--threads', type=int, default=8)
    p.add_argument('--samtools', default='samtools'); p.set_defaults(f=cmd_count)
    p = sub.add_parser('call'); p.add_argument('samples'); p.add_argument('gtf'); p.add_argument('target'); p.add_argument('outdir'); p.add_argument('--de'); p.add_argument('--strict', action='store_true'); p.set_defaults(f=cmd_call)
    a = ap.parse_args(); a.f(a)


if __name__ == '__main__':
    main()
