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

mkdir ${OUTDIR}/2_SVs/Delly/${SAMPLE}_Results/

# default deletions
delly call \
-q 30 \
-g ${REF} \
-o ${OUTDIR}/2_SVs/Delly/${SAMPLE}_Results/${SAMPLE}.DEL.bcf \
${OUTDIR}/1_Mapping/${SAMPLE}.RG.bam

# insertions
delly call \
-q 30 \
-t INS \
-g ${REF} \
-o ${OUTDIR}/2_SVs/Delly/${SAMPLE}_Results/${SAMPLE}.INS.bcf \
${OUTDIR}/1_Mapping/${SAMPLE}.RG.bam

# duplications
delly call \
-q 30 \
-t DUP \
-g ${REF} \
-o ${OUTDIR}/2_SVs/Delly/${SAMPLE}_Results/${SAMPLE}.DUP.bcf \
${OUTDIR}/1_Mapping/${SAMPLE}.RG.bam

# merge Del and Ins
bcftools merge \
-o ${OUTDIR}/2_SVs/Delly/${SAMPLE}_Results/${SAMPLE}.bcf \
-f PASS \
-O b \
--force-samples \
${OUTDIR}/2_SVs/Delly/${SAMPLE}_Results/${SAMPLE}.DEL.bcf \
${OUTDIR}/2_SVs/Delly/${SAMPLE}_Results/${SAMPLE}.INS.bcf
