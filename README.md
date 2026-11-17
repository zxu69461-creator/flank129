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

Ready-made lists are in `sites/`, so this step can be skipped:

- `sites_129_shared.GRCm38.tsv.gz`: the 4,458,459 sites above (GRCm38, as used in the paper).
- `sites_129_shared.GRCm39.tsv.gz`: 4,456,540 of them in GRCm39 coordinates. They were converted with the UCSC mm10-to-mm39 chain and
  kept only where the GRCm39 reference base is the site's C57BL/6 base. Sites in segments whose orientation differs between the
  assemblies (mostly on chromosome 14) are given with their alleles on the GRCm39 forward strand. `sites_GRCm39_report.tsv` counts the
  1,919 sites dropped, and `tools/convert_sites_grcm39.py` is the conversion.

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

Rules: sites with at least 4 reads; a 1 Mb bin with at least 3 such sites is 129-carrying when its 129 allele fraction is at least 0.2;
a bin enters the comparison of groups only when every sample of both groups has a call there; adjacent separating bins are merged;
segments within 20 Mb of the target (or containing it) are reported as near the target.

Reading the result: a knockout-only 129 segment near the target means that genes inside it differ between the groups in strain as well
as in genotype at the target; their expression changes cannot be attributed to the knockout from these data alone. Controls that carry
the same flank (littermates heterozygous for the allele, or floxed controls in conditional designs) remove the difference.

## Example

`example/samples_GSE173144.tsv` is the samples table of the *Cxcr5* study shown in the paper (three knockout and three wild-type mice;
the count files come from step 2). `example/segments_GSE173144.tsv` is the output of `call` for it: knockout-only 129 DNA in eight
segments between 32 and 76 Mb of chromosome 9, around *Cxcr5*.

## Resources

Calling a six-sample study from its count files takes about 2 s and 0.5 GB of memory on a laptop. The cost lies in aligning the reads:
in the paper, downloading, aligning and counting took about 16 machine-minutes per run on eight-core cloud machines.

## Licence and citation

MIT licence (see `LICENSE`). If you use flank129, please cite the article in `CITATION.cff`.
