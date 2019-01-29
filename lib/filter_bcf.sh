#!/bin/bash
#==============================================================================#
# Filters the results from each of the variant calling tools, saving calls that pass and
# fail into separate bcf files
# USAGE:
# -s [samples_names.txt] list of sample names
# -o [path/tp/outdir] output working directory for Varsect project
# [-D] [-F] [-D] [-P] Variant calling tool results to filter
#==============================================================================#

while getopts "s:r:o:DFGP" flag; do
  case "${flag}" in
    s ) SAMPLE="${OPTARG}" ;;
    o ) OUTDIR="${OPTARG}" ;;
    D ) DELLY=true ;;
    F ) FREEBAYES=true ;;
    G ) GATK=true ;;
    P ) PINDEL=true ;;
  esac
done


# Delly
if [ "${DELLY}" = "true" ] ; then
    echo "Filtering Delly results for ${SAMPLE}"

    # merge Del and Ins (has also been put in delly.sh)
    bcftools merge \
    -o ${OUTDIR}/2_SVs/Delly/${SAMPLE}_Results/${SAMPLE}.bcf \
    -f PASS \
    -O b \
    --force-samples \
    ${OUTDIR}/2_SVs/Delly/${SAMPLE}_Results/${SAMPLE}.DEL.bcf \
    ${OUTDIR}/2_SVs/Delly/${SAMPLE}_Results/${SAMPLE}.INS.bcf

    # filter for min support and remove LowQual
    bcftools view \
    -i 'MIN( ((FMT/RV) / ((FMT/RR)+(FMT/RV))) > 0.7)' \
    ${OUTDIR}/2_SVs/Delly/${SAMPLE}_Results/${SAMPLE}.bcf | \
    grep -v LowQual | \
    bcftools view \
    -O b \
    > ${OUTDIR}/2_SVs/Delly/${SAMPLE}_Results/${SAMPLE}.PASS.bcf

    # failed (low support)
    bcftools view \
    -i 'MIN( ((FMT/RV) / ((FMT/RR)+(FMT/RV))) < 0.7)' \
    ${OUTDIR}/2_SVs/Delly/${SAMPLE}_Results/${SAMPLE}.bcf \
    > ${OUTDIR}/2_SVs/Delly/${SAMPLE}_Results/tmp.vcf

    # failed (LowQual)
    bcftools view \
    -i 'MIN( ((FMT/RV) / ((FMT/RR)+(FMT/RV))) > 0.7)' \
    ${OUTDIR}/2_SVs/Delly/${SAMPLE}_Results/${SAMPLE}.bcf | \
    grep LowQual | \
    grep -v '^#' \
    >> ${OUTDIR}/2_SVs/Delly/${SAMPLE}_Results/tmp.vcf

    # convert tmp to .bcf format
    bcftools view \
    -o ${OUTDIR}/2_SVs/Delly/${SAMPLE}_Results/${SAMPLE}.FAIL.bcf \
    -O b \
    ${OUTDIR}/2_SVs/Delly/${SAMPLE}_Results/tmp.vcf

    rm ${OUTDIR}/2_SVs/Delly/${SAMPLE}_Results/tmp.vcf
fi

# Freebayes
if [ "${FREEBAYES}" = "true" ] ; then
    echo "Filtering Freebayes results for ${SAMPLE}"

    bcftools view \
    -o ${OUTDIR}/2_SVs/Freebayes/${SAMPLE}_Results/${SAMPLE}.PASS.bcf \
    -O b \
    -i 'MIN( ((FMT/AO[0]) / ((FMT/AO[0])+(FMT/RO[0]))) > 0.7)' \
    ${OUTDIR}/2_SVs/Freebayes/${SAMPLE}_Results/${SAMPLE}.bcf

    bcftools view \
    -o ${OUTDIR}/2_SVs/Freebayes/${SAMPLE}_Results/${SAMPLE}.FAIL.bcf \
    -O b \
    -i 'MIN( ((FMT/AO[0]) / ((FMT/AO[0])+(FMT/RO[0]))) < 0.7)' \
    ${OUTDIR}/2_SVs/Freebayes/${SAMPLE}_Results/${SAMPLE}.bcf
fi


# GATK
if [ "${GATK}" = "true" ] ; then
    echo "Filtering GATK results for ${SAMPLE}"

    bcftools view \
    -o ${OUTDIR}/2_SVs/GATK/${SAMPLE}_Results/tmp.vcf \
    ${OUTDIR}/2_SVs/GATK/${SAMPLE}_Results/${SAMPLE}.bcf

    gatk-launch IndexFeatureFile --feature-file ${OUTDIR}/2_SVs/GATK/${SAMPLE}_Results/tmp.vcf

    gatk-launch VariantFiltration \
    --output ${OUTDIR}/2_SVs/GATK/${SAMPLE}_Results/tmp.filtered.vcf \
    --variant ${OUTDIR}/2_SVs/GATK/${SAMPLE}_Results/tmp.vcf \
    --cluster-size 3 --cluster-window-size 10 \
    --filter-expression "QD < 10.0" --filter-name "FAIL-QD" \
    --filter-expression "MQ < 30.0" --filter-name "FAIL-MQ" \
    --filter-expression "FS > 10.0" --filter-name "FAIL-FS" \
    --filter-expression "QUAL < 30.0" --filter-name "FAIL-QUAL" \
    --genotype-filter-expression "DP < 10" --genotype-filter-name "FAIL-DP" \
    --genotype-filter-expression "AD.1 < (DP*0.7)" --genotype-filter-name "FAIL-low_freq"


    # put PASS and FAILED calls in separate files
    bcftools view \
    -h \
    ${OUTDIR}/2_SVs/GATK/${SAMPLE}_Results/tmp.filtered.vcf \
    > ${OUTDIR}/2_SVs/GATK/${SAMPLE}_Results/${SAMPLE}.PASS.vcf

    bcftools view \
    -H \
    ${OUTDIR}/2_SVs/GATK/${SAMPLE}_Results/tmp.filtered.vcf | \
    grep -v FAIL >> ${OUTDIR}/2_SVs/GATK/${SAMPLE}_Results/${SAMPLE}.PASS.vcf

    bcftools view \
    -h \
    ${OUTDIR}/2_SVs/GATK/${SAMPLE}_Results/tmp.filtered.vcf \
    > ${OUTDIR}/2_SVs/GATK/${SAMPLE}_Results/${SAMPLE}.FAIL.vcf

    bcftools view \
    -H \
    ${OUTDIR}/2_SVs/GATK/${SAMPLE}_Results/tmp.filtered.vcf | \
    grep FAIL >> ${OUTDIR}/2_SVs/GATK/${SAMPLE}_Results/${SAMPLE}.FAIL.vcf

    bcftools view \
    -o ${OUTDIR}/2_SVs/GATK/${SAMPLE}_Results/${SAMPLE}.PASS.bcf \
    -O b \
    ${OUTDIR}/2_SVs/GATK/${SAMPLE}_Results/${SAMPLE}.PASS.vcf

    bcftools view \
    -o ${OUTDIR}/2_SVs/GATK/${SAMPLE}_Results/${SAMPLE}.FAIL.bcf \
    -O b \
    ${OUTDIR}/2_SVs/GATK/${SAMPLE}_Results/${SAMPLE}.FAIL.vcf

    rm ${OUTDIR}/2_SVs/GATK/${SAMPLE}_Results/*.vcf*
fi

# Pindel
if [ "${PINDEL}" = "true" ] ; then
    echo "Filtering Pindel results for ${SAMPLE}"

    bcftools view \
    -o ${OUTDIR}/2_SVs/Pindel/${SAMPLE}_Results/${SAMPLE}.PASS.bcf \
    -O b \
    -i 'MIN( ((FMT/AD[1]) / ((FMT/AD[1])+(FMT/AD[0]))) > 0.7)' \
    ${OUTDIR}/2_SVs/Pindel/${SAMPLE}_Results/${SAMPLE}.bcf

    bcftools view \
    -o ${OUTDIR}/2_SVs/Pindel/${SAMPLE}_Results/${SAMPLE}.FAIL.bcf \
    -O b \
    -i 'MIN( ((FMT/AD[1]) / ((FMT/AD[1])+(FMT/AD[0]))) < 0.7) & MIN( (FMT/AD[1]) / ((FMT/AD[1])+(FMT/AD[0])) > 0.2)' \
    ${OUTDIR}/2_SVs/Pindel/${SAMPLE}_Results/${SAMPLE}.bcf
fi
