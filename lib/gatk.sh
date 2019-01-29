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

mkdir ${OUTDIR}/2_SVs/GATK/${SAMPLE}_Results/

gatk-launch HaplotypeCaller \
-I ${OUTDIR}/1_Mapping/${SAMPLE}.RG.bam \
-O ${OUTDIR}/2_SVs/GATK/${SAMPLE}_Results/${SAMPLE}.vcf \
-R ${REF} \
-ploidy 1 \
--read-filter MappingQualityReadFilter \
--minimum-mapping-quality 30

bcftools view \
-o ${OUTDIR}/2_SVs/GATK/${SAMPLE}_Results/${SAMPLE}.bcf \
-O b \
${OUTDIR}/2_SVs/GATK/${SAMPLE}_Results/${SAMPLE}.vcf

rm ${OUTDIR}/2_SVs/GATK/${SAMPLE}_Results/${SAMPLE}.vcf*
