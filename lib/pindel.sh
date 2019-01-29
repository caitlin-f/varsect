#!/bin/bash
# t = threads
# s = list of sample names
# r = reference fasta file
# o = output directory

# Pindel Arguments
# -M minimum support for event (default 1)
# -A minimum anchor quality (default 0)
# -l report long insertions (default false)
# -x change maximum size of structural variant to detect to 3 = 2048bp (default 2 = 512bp)
# -r report inversions (default true), -v min inversion size (default 50)
# -q detect dispersed duplications (default false)

while getopts "t:s:r:1:2:o:" flag; do
  case "${flag}" in
    t ) T="${OPTARG}" ;;
    s ) SAMPLE="${OPTARG}" ;;
    r ) REF="${OPTARG}" ;;
    o ) OUTDIR="${OPTARG}" ;;
  esac
done

mkdir -p ${OUTDIR}/2_SVs/Pindel/${SAMPLE}_Results/tmp/

picard CollectInsertSizeMetrics \
I=${OUTDIR}/1_Mapping/${SAMPLE}.RG.bam \
O=${OUTDIR}/1_Mapping/${SAMPLE}.picard.metrics \
H=${OUTDIR}/1_Mapping/${SAMPLE}.picard.hist

INSERT=$(head -n 8 ${OUTDIR}/1_Mapping/${SAMPLE}.picard.metrics | tail -n 1 | cut -f 6)
# make config file
printf "${OUTDIR}/1_Mapping/${SAMPLE}.RG.bam\t${INSERT}\t${SAMPLE}" > ${OUTDIR}/2_SVs/Pindel/${SAMPLE}_Results/${SAMPLE}.pd.conf

pindel \
-T ${T} \
-l \
-M 4 \
-A 30 \
-x 3 \
-f ${REF} \
-i ${OUTDIR}/2_SVs/Pindel/${SAMPLE}_Results/${SAMPLE}.pd.conf \
-o ${OUTDIR}/2_SVs/Pindel/${SAMPLE}_Results/tmp/${SAMPLE}.pd.out

# Convert files to pindel output to bcf format
for SV in D SI LI TD
do
    pindel2vcf \
    -p ${OUTDIR}/2_SVs/Pindel/${SAMPLE}_Results/tmp/${SAMPLE}.pd.out_${SV} \
    -r ${REF} \
    -R ${REF} \
    -d date \
    -v ${OUTDIR}/2_SVs/Pindel/${SAMPLE}_Results/tmp/${SAMPLE}.${SV}.vcf

    # Pindel2vcf output does not convert to bcf simply
    # need to perform following to compress to bcf
    bgzip --force ${OUTDIR}/2_SVs/Pindel/${SAMPLE}_Results/tmp/${SAMPLE}.${SV}.vcf
    bcftools index -f ${OUTDIR}/2_SVs/Pindel/${SAMPLE}_Results/tmp/${SAMPLE}.${SV}.vcf.gz

    # Can now convert to bcf
    bcftools view \
    -o ${OUTDIR}/2_SVs/Pindel/${SAMPLE}_Results/${SAMPLE}.${SV}.bcf \
    -O b \
    ${OUTDIR}/2_SVs/Pindel/${SAMPLE}_Results/tmp/${SAMPLE}.${SV}.vcf.gz

    bcftools index -f ${OUTDIR}/2_SVs/Pindel/${SAMPLE}_Results/${SAMPLE}.${SV}.bcf
done

bcftools merge \
-o ${OUTDIR}/2_SVs/Pindel/${SAMPLE}_Results/${SAMPLE}.bcf \
-f PASS \
-O b \
--force-samples \
${OUTDIR}/2_SVs/Pindel/${SAMPLE}_Results/${SAMPLE}.D.bcf \
${OUTDIR}/2_SVs/Pindel/${SAMPLE}_Results/${SAMPLE}.SI.bcf

# remove intermediate directory
rm -r ${OUTDIR}/2_SVs/Pindel/${SAMPLE}_Results/tmp
