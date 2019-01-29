#!/bin/bash
# s = list of sample names
# r = reference fasta file
# o = output directory

while getopts "s:r:o:" flag; do
  case "${flag}" in
    s ) SAMPLE="${OPTARG}" ;;
    r ) REF="${OPTARG}" ;;
    o ) OUTDIR="${OPTARG}" ;;
  esac
done

mkdir ${OUTDIR}/2_SVs/Freebayes/${SAMPLE}_Results/

freebayes -f ${REF} \
-p 1 \
-C 4 \
-m 30 \
${OUTDIR}/1_Mapping/${SAMPLE}.RG.bam | \
bcftools view \
-o ${OUTDIR}/2_SVs/Freebayes/${SAMPLE}_Results/${SAMPLE}.bcf \
-O b
