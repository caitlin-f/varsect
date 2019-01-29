# Varsect

![Varsect logo](lib/varsect_logo_2.png)

Varsect (VARiant calling by interSECTion) is a pipeline designed to report single nucleotide polymorphisms (SNPs),
small insertions and deletions (indels) and large variants (specifically large
deletions and translocations) in microbial genomes.

Varsect has been optimised for use with the Slurm workflow manager in
Linux environments to efficiently analyse large datasets.

Varsect uses existing freely available bioinformatics tools for mapping and variant calling including BWA, Delly, GATK, Picard, Pindel and Samtools. Results of the various variant calling tools is collated and only high confidence and high frequency variants are reported.

Output files include vcf files suitable for variant annotation, sample x variant csv files, and alignment files (both phylip and nexus) for phylogenomic analysis.

### Usage:
varsect_batch.py [-h] -r reference.fa -o path/to/outdir -s samples.txt -t threads [-E recomb.gff] -A
