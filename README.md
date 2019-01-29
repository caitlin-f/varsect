# Varsect

![Varsect logo](https://github.com/caitlin-f/varsect/lib/varsect_logo.png)

Varsect is a pipeline designed to report single nucleotide polymorphisms (SNPs),
small insertions and deletions (indels) and large variants, specifically large
deletions and translocations in microbial genomes.

Varsect has been optimised for use with the Slurm workflow manager in
Linux environments to efficiently analyse large datasets.

Varsect uses existing freely available bioinformatics tools for mapping and
variant calling including BWA, Delly, GATK, Picard, Pindel and Samtools.

Usage: varsect_batch.py [-h] -r [reference.fa] -o [path/to/outdir] -s  
                        [samples.txt] -t [threads] [-E [recomb.gff]] [-A] [-M]  
                        [-D] [-F] [-G] [-P] [-C] [-I] [-R] [-B]  

optional arguments: \n
  -h, --help           show this help message and exit  
  -r [reference.fa]    Reference file in .fasta format  
  -o [path/to/outdir]  Output directory  
  -s [samples.txt]     Read filenames to analyse  
  -t [threads]         Number of threads  
  -E [recomb.gff]      Remove recombination regions specified in a gff file  
  -A                   Run all steps  
  -M                   Perform Mapping  
  -D                   Run Delly  
  -F                   Run Freebayes  
  -G                   Run GATK-HaplotypeCaller  
  -P                   Run Pindel  
  -C                   Collate data  
  -I                   Make Samplot Images  
  -R                   Run RaxML  
  -B                   Run MrBayes  
