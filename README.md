# Varsect

![Varsect logo](lib/varsect_logo_2.png)

Varsect *(VARiant calling by interSECTion)* is a pipeline designed to report single nucleotide polymorphisms (SNPs),
small insertions and deletions (indels) and large variants (specifically large
deletions and translocations) in microbial genomes.

Varsect has been optimised for use with the Slurm workflow manager in
Linux environments to efficiently analyse large datasets.

Varsect uses existing freely available bioinformatics tools for mapping and variant calling including BWA, Delly, GATK, Picard, Pindel and Samtools. Results of the various variant calling tools is collated and only high confidence and high frequency variants are reported.

Output files include vcf files suitable for variant annotation, sample x variant csv files, and alignment files (both phylip and nexus) for phylogenomic analysis.

### Usage:
```
varsect_batch.py [-h] -r reference.fa -o path/to/outdir -s samples.txt -t threads [-E recomb.gff] [-A]


Arguments:
  -h, --help       show this help message and exit
  -r STR           reference file in fasta format (required)
  -o STR           output directory (required)
  -s FILE          line separated list of paired-end fastq.gz filenames (required)
  -t INT           number of threads (required)
  -E FILE          remove recombination regions specified in a gff file
  -A               write batch scripts for all stages and all tools

The following may be specified instead of -A to write batch scripts for individual tools or steps:
  -M               mapping with BWA mem
  -D               call large variants with Delly
  -F               call SNPs and small indels with Freebayes
  -G               calls SNPs and small indels with GATK-HaplotypeCaller
  -P               call large variants with Pindel
  -C               collate data (filtering, intersection, collation)
  -I               Samplot images
  -R               RaxML
  -B               MrBayes
```
