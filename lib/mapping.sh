#!/bin/bash
# t = threads
# s = list of sample names
# r = reference fasta file
# 1 = list of first in pair fq/fq.qz reads
# 2 = list of second in pair fq/fq.gz reads
# o = output directory

while getopts ":t:s:r:1:2:o:" flag; do
  case "${flag}" in
    t ) T="${OPTARG}" ;;
    s ) SAMPLE="${OPTARG}" ;;
    r ) REF="${OPTARG}" ;;
    1 ) FWD="${OPTARG}" ;;
    2 ) REV="${OPTARG}" ;;
    o ) OUTDIR="${OPTARG}" ;;
  esac
done


bwa mem -t ${T} ${REF} ${FWD} ${REV} | \
samtools sort | \
samtools view -b > ${OUTDIR}/1_Mapping/${SAMPLE}.sorted.bam

picard MarkDuplicates \
I=${OUTDIR}/1_Mapping/${SAMPLE}.sorted.bam \
O=${OUTDIR}/1_Mapping/${SAMPLE}.noRG.bam \
M=${OUTDIR}/1_Mapping/${SAMPLE}.markdups.metrics \

rm ${OUTDIR}/1_Mapping/${SAMPLE}.sorted.bam*

samtools index ${OUTDIR}/1_Mapping/${SAMPLE}.noRG.bam

picard AddOrReplaceReadGroups \
I=${OUTDIR}/1_Mapping/${SAMPLE}.noRG.bam \
O=${OUTDIR}/1_Mapping/${SAMPLE}.RG.bam \
RGID=${SAMPLE} RGLB=lib${SAMPLE} RGPL=illumina RGPU=unit${SAMPLE} RGSM=SM${SAMPLE}

rm ${OUTDIR}/1_Mapping/${SAMPLE}.noRG.bam*

samtools index ${OUTDIR}/1_Mapping/${SAMPLE}.RG.bam

samtools mpileup -aa -q 30 ${OUTDIR}/1_Mapping/${SAMPLE}.RG.bam | \
awk '($4 < 11)' > ${OUTDIR}/1_Mapping/${SAMPLE}_low_cov.mpileup
