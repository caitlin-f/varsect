# Varsect

![Varsect logo](lib/varsect_logo.png)

Varsect is a pipeline designed to report single nucleotide polymorphisms (SNPs),
small insertions and deletions (indels) and large variants, specifically large
deletions and translocations in microbial genomes.

Varsect has been optimised for use with the Slurm workflow manager in
Linux environments to efficiently analyse large datasets.

Varsect uses existing freely available bioinformatics tools for mapping and
variant calling including BWA, Delly, GATK, Picard, Pindel and Samtools.

Usage: varsect_batch.py [-h] -r [reference.fa] -o [path/to/outdir] -s [samples.txt] -t [threads] [-E [recomb.gff] [-A] 

optional arguments:  
&nbsp;&nbsp;-h, --help&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;show this help message and exit  
&nbsp;&nbsp;-r [reference.fa]&nbsp;&nbsp;&nbsp;&nbsp;Reference file in .fasta format  
&nbsp;&nbsp;-o [path/to/outdir]&nbsp;&nbsp;Output directory  
&nbsp;&nbsp;-s [samples.txt]&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Read filenames to analyse  
&nbsp;&nbsp;-t [threads]&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Number of threads  
&nbsp;&nbsp;-E [recomb.gff]&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Remove recombination regions specified in a gff file  
&nbsp;&nbsp;-A&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Run all steps  
&nbsp;&nbsp;-M&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Perform Mapping  
&nbsp;&nbsp;-D&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Run Delly  
&nbsp;&nbsp;-F&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Run Freebayes  
&nbsp;&nbsp;-G&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Run GATK-HaplotypeCaller  
&nbsp;&nbsp;-P&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Run Pindel  
&nbsp;&nbsp;-C&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Collate data  
&nbsp;&nbsp;-I&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Make Samplot Images  
&nbsp;&nbsp;-R&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Run RaxML  
&nbsp;&nbsp;-B&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Run MrBayes  
