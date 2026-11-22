# flank129

Reads the 129-strain DNA around a targeted gene from a knockout study's own RNA-seq reads, and reports whether the knockouts carry 129
DNA there that their controls lack.

Most knockout alleles were made in 129 embryonic stem cells and backcrossed to C57BL/6. Knockouts then differ from their controls not
only in the targeted gene but in the 129 DNA around it, and genes in that flank can change expression because of the strain, not the
knockout. flank129 makes this visible for any study with raw reads.

Requirements: Python 3 with pandas and numpy; samtools (for `count`); bcftools (for `sites`); an aligner such as HISAT2 with a GRCm38
or GRCm39 index.

## 1. Sites (once)

The sites are the single-nucleotide positions where 129S1/SvImJ, 129P2/OlaHsd and 129S5/SvEvBrd are all homozygous for the same
non-reference allele (Mouse Genomes Project release REL-1505, GRCm38; 4,458,459 sites):

    flank129.py sites 129S1_SvImJ.mgp.v5.snps.dbSNP142.vcf.gz 129P2_OlaHsd.mgp.v5.snps.dbSNP142.vcf.gz \
        129S5SvEvBrd.mgp.v5.snps.dbSNP142.vcf.gz sites_129_shared.tsv.gz

The VCFs are at https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-1505-SNPs_Indels/strain_specific_vcfs/.

Other donor strains: give that strain's VCF alone, e.g. for alleles edited in (C57BL/6 x DBA/2) zygotes

    flank129.py sites DBA_2J.mgp.v5.snps.dbSNP142.vcf.gz sites_DBA_2J.tsv.gz

and use that list in steps 2 and 3. The 129 sites are not specific to 129: at about half of them other strains carry the same allele,
so on a C57BL/6 background they mark non-C57BL/6 DNA in general.

Ready-made lists are in `sites/`, so this step can be skipped:

- `sites_129_shared.GRCm38.tsv.gz`: the 4,458,459 sites above (GRCm38, as used in the paper).
- `sites_129_shared.GRCm39.tsv.gz`: 4,456,540 of them in GRCm39 coordinates. They were converted with the UCSC mm10-to-mm39 chain and
  kept only where the GRCm39 reference base is the site's C57BL/6 base. Sites in segments whose orientation differs between the
  assemblies (mostly on chromosome 14) are given with their alleles on the GRCm39 forward strand. `sites_GRCm39_report.tsv` counts the
  1,919 sites dropped, and `tools/convert_sites_grcm39.py` is the conversion.
- `sites_CBA_J.GRCm38.tsv.gz`, `sites_DBA_2J.GRCm38.tsv.gz`, `sites_BALB_cJ.GRCm38.tsv.gz`, `sites_C3H_HeJ.GRCm38.tsv.gz`: the
  homozygous non-reference single-nucleotide sites (FILTER PASS) of each strain in the same release, for alleles edited in hybrid zygotes
  or made in stem cells of these strains (`tools/convert_sites_grcm39.py` converts them to GRCm39).
- `sites_<strain>_private.GRCm38.tsv.gz` (129_shared, CBA_J, DBA_2J, BALB_cJ, C3H_HeJ): the subset of each list where no other of these
  five carries the alternative allele (1,364,396 sites for 129; 231,937–935,659 for the others). Many alleles are shared between
  laboratory strains, so a CBA/J segment also lights up the 129 list and the reverse; calling a study with two private lists tells the
  donor strain apart from 129 (see Validation).

Use the list that matches the assembly your reads are aligned to, with a GTF of the same assembly in step 3.

## 2. Counts per run

Align every run to GRCm38 (or GRCm39, with the GRCm39 list) and count the reads with the C57BL/6 and the 129 base at each site:

    hisat2 -p 8 --no-unal -x grcm38/genome -U run.fastq | samtools sort -o run.bam - && samtools index run.bam
    flank129.py count run.bam sites_129_shared.tsv.gz run.ad.tsv.gz

## 3. Segments for one study

`samples.tsv` lists every file with its sample and group (several runs of one sample are summed):

    sample      group     file
    GSM5261098  knockout  GSM5261098_SRR14304915.ad.tsv.gz
    GSM5261095  control   GSM5261095_SRR14304912.ad.tsv.gz

    flank129.py call example/samples_GSE173144.tsv Mus_musculus.GRCm38.102.gtf.gz Cxcr5 out/ [--de de.tsv]

Output: `out/bins.tsv.gz` (every sample and 1 Mb bin: sites, reads, 129 allele fraction, call) and `out/segments.tsv` (segments in which
every knockout and no control sample carries 129 alleles, and the reverse, with their distance to the target). With `--de` (columns
`gene`, `fdr`), the genes differentially expressed inside the knockout-only segments near the target are counted against all others.

Rules (version 1.1): sites with at least 4 reads; a 1 Mb bin with at least 3 such sites is 129-carrying when its 129 allele fraction
is at least 0.2. A bin enters the comparison of groups when each group has calls in at least max(2, half) of its samples, and separates
the groups when every called knockout and no called control is 129-carrying (or the reverse). Separating bins are joined across gaps of
up to 5 bins that lack enough reads, but not across a bin that was compared and did not separate. Segments within 20 Mb of the target
(or containing it) are reported as near the target; `bridged_mb` gives the length of the joined gaps. `--strict` applies the rule of
version 1.0 (every sample called, only adjacent bins merged), which split long segments wherever one sample lacked reads.

Reading the result: a knockout-only 129 segment near the target means that genes inside it differ between the groups in strain as well
as in genotype at the target; their expression changes cannot be attributed to the knockout from these data alone. Controls that carry
the same flank (littermates heterozygous for the allele, or floxed controls in conditional designs) remove the difference.

## Validation

Bin calls from RNA-seq reads were checked against DNA genotypes of the same strains, processed exactly as above:

- 129 sites: kidney RNA-seq of 16 Collaborative Cross strains (one male and one female each; GEO GSE199702) against founder mosaics
  from whole-genome sequencing, at bins where the animals' MiniMUGA genotypes agree. Bins called 129-carrying: sensitivity 99.7% in
  129S1/129S1 bins, specificity 98.9% in C57BL/6J/C57BL/6J bins (99.8% and 99.2% after thinning to 5 million reads, which calls
  about 850 instead of 1,200 bins per sample). Bins of the other founders, which share many alleles with 129, are called 129-carrying
  in about two thirds of cases, so a segment called from the 129 list means non-C57BL/6 DNA.
- DBA/2J sites: hypothalamus RNA-seq of 16 BXD strains (GSE79551) against their marker genotypes: sensitivity 98.8%, specificity 96.2%.
- Segment ends fall a median of 0.5 Mb from the DNA boundaries (90th percentile 7 Mb).
- Inbred controls of the same study (C57BL/6J, CBA/J, DBA/2J, BALB/cJ, C3H/HeJ): C57BL/6J is called carrying in at most 0.3% of bins
  at any strain's list. At sites private to one strain among 129, CBA/J, DBA/2J, BALB/cJ and C3H/HeJ, each strain lights up its own
  private sites (92–99% of bins) and almost none of the others (at most 1%), which separates a donor strain from 129.

## Example

`example/samples_GSE173144.tsv` is the samples table of the *Cxcr5* study shown in the paper (three knockout and three wild-type mice;
the count files come from step 2). `example/segments_GSE173144.tsv` is the output of `call` for it: one knockout-only 129 segment of
44 Mb, from 32 to 76 Mb of chromosome 9, containing *Cxcr5* (version 1.0 reported it as eight pieces), and a 1 Mb segment on
chromosome 11.

## Resources

Calling a six-sample study from its count files takes about 2 s and 0.5 GB of memory on a laptop. The cost lies in aligning the reads:
in the paper, downloading, aligning and counting took about 16 machine-minutes per run on eight-core cloud machines.

## Licence and citation

Archived at Zenodo: https://doi.org/10.5281/zenodo.24309333 (all versions; 1.0.0 is https://doi.org/10.5281/zenodo.24309334). MIT licence (see `LICENSE`). If you use flank129, please cite the article in `CITATION.cff`.
