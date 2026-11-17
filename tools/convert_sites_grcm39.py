"""GRCm39 version of the 129 site list used by flank129 (2026-11-17, for users who align to GRCm39).

The 4,458,459 sites where 129S1/SvImJ, 129P2/OlaHsd and 129S5/SvEvBrd share a homozygous non-reference allele (Mouse Genomes Project
REL-1505, GRCm38) are converted with the UCSC mm10-to-mm39 chain. A site is kept when it converts to a single position on the same
chromosome and the GRCm39 reference base there is the site's reference base. Sites in segments whose orientation differs between the
assemblies (most of them in one part of chromosome 14) map to the minus strand; their alleles are complemented to the GRCm39 forward
strand. Sites whose GRCm39 base is the 129 allele (the reference changed) or anything else are dropped and counted.

Usage: 23_sites_grcm39.py SITES_GRCm38.tsv.gz CHAIN.over.chain.gz GRCm39.primary_assembly.fa.gz OUT_DIR
Writes OUT_DIR/sites_129_shared.GRCm38.tsv.gz (copy), OUT_DIR/sites_129_shared.GRCm39.tsv.gz and OUT_DIR/sites_GRCm39_report.tsv."""
import gzip
import shutil
import sys
from pathlib import Path
import pandas as pd
from pyliftover import LiftOver

src, chain, fasta, out = sys.argv[1], sys.argv[2], sys.argv[3], Path(sys.argv[4]); out.mkdir(parents=True, exist_ok=True)
S = pd.read_csv(src, sep='\t', dtype={'chrom': str})
lo = LiftOver(chain)
COMP = str.maketrans('ACGT', 'TGCA')
new_pos, status, flip = [], [], []
for c, p in zip(S.chrom.values, S.pos.values):
    r = lo.convert_coordinate('chr' + c, int(p) - 1)          # 0-based in, 0-based out (pyliftover gives the forward-strand position)
    if not r: new_pos.append(-1); status.append('unmapped'); flip.append(False); continue
    if len(r) > 1: new_pos.append(-1); status.append('multiple'); flip.append(False); continue
    ch, q, strand, _ = r[0]
    if ch != 'chr' + c: new_pos.append(-1); status.append('other chromosome'); flip.append(False); continue
    new_pos.append(q + 1); status.append('lifted'); flip.append(strand == '-')
S['pos39'] = new_pos; S['status'] = status; S['flip'] = flip
S.loc[S.flip, 'ref'] = S.loc[S.flip, 'ref'].str.translate(COMP)        # alleles on the GRCm39 forward strand
S.loc[S.flip, 'alt'] = S.loc[S.flip, 'alt'].str.translate(COMP)
print('minus-strand sites (alleles complemented):', int(S.flip.sum()), flush=True)
print(S.status.value_counts().to_dict(), flush=True)

# reference base check, one chromosome at a time from the gzipped FASTA
want = {c: g for c, g in S[S.status == 'lifted'].groupby('chrom')}
base = pd.Series('', index=S.index)


def check(name, seq):
    g = want.get(name)
    if g is None: return
    s = ''.join(seq).upper()
    base.loc[g.index] = [s[p - 1] if 0 < p <= len(s) else '' for p in g.pos39.values]
    print(f'chr{name}: {len(g):,} sites checked', flush=True)


name, seq = None, []
with gzip.open(fasta, 'rt') as f:
    for line in f:
        if line.startswith('>'):
            if name is not None: check(name, seq)
            name, seq = line[1:].split()[0], []
        else:
            seq.append(line.rstrip())
    if name is not None: check(name, seq)
S['ref39'] = base
lifted = S.status == 'lifted'
S.loc[lifted & (S.ref39 == S.alt), 'status'] = 'reference is the 129 allele'
S.loc[lifted & (S.ref39 != S.ref) & (S.ref39 != S.alt), 'status'] = 'other base'
S.loc[lifted & (S.ref39 == S.ref), 'status'] = 'kept'
rep = S.status.value_counts().rename_axis('status').reset_index(name='sites'); rep.to_csv(out / 'sites_GRCm39_report.tsv', sep='\t', index=False)
print(rep.to_string(index=False), flush=True)
K = S[S.status == 'kept'][['chrom', 'pos39', 'ref', 'alt']].rename(columns={'pos39': 'pos'}).sort_values(['chrom', 'pos'])
K.to_csv(out / 'sites_129_shared.GRCm39.tsv.gz', sep='\t', index=False)
shutil.copyfile(src, out / 'sites_129_shared.GRCm38.tsv.gz')
print(f'{len(K):,} GRCm39 sites -> {out}')
