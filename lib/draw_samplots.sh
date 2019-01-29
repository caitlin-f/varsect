#!/bin/bash

while getopts "n:o:c:t:s:e:" flag; do
  case "${flag}" in
    n ) SAMPLE="${OPTARG}" ;;
    o ) OUTDIR="${OPTARG}" ;;
    c ) CHROM="${OPTARG}" ;;
    t ) SV="${OPTARG}" ;;
    s ) POS="${OPTARG}" ;;
    e ) END="${OPTARG}" ;;
  esac
done

samplot.py \
-n ${SAMPLE} \
-b ${OUTDIR}/1_Mapping/${SAMPLE}.RG.bam \
-o ${OUTDIR}/4_Images/${SAMPLE}_Results/${SAMPLE}.${SV}.${POS}.png \
-c ${CHROM} \
-t ${SV} \
-s ${POS} \
-e ${END}

IFS=$'\n'
for result in $( bcftools query -f "%CHROM %POS %INFO/END %INFO/SVTYPE\n" -i 'INFO/SVTYPE = "DEL"' ${DIR}/2_SVs/unique_pos.vcf ); do
  IFS=$' '
  sv=($result)
  # echo CHROM ${sv[0]}
  # echo POS ${sv[1]}
  # echo END ${sv[2]}
  # echo TYPE ${sv[3]}
  samplot.py \
      -n ${samples} \
      -b ${bam_files} \
      -o ${DIR}/4_Images/${sv[3]}_${sv[1]}.png \
      -c ${sv[0]} \
      -s ${sv[1]} \
      -e ${sv[2]} \
      -t ${sv[3]} \
      -d 100
  echo Done ${sv[0]} ${sv[1]} ${sv[3]}
  IFS=$'\n'
done

for result in $( bcftools query -f "%CHROM %POS %INFO/END %INFO/SVTYPE %INFO/LEN \n" -i 'INFO/SVTYPE = "TSLN"' ${DIR}/2_SVs/unique_pos.vcf ); do
  IFS=$' '
  sv=($result)
  e1=$(( ${sv[1]}+1 ))
  samplot.py \
      -n ${samples} \
      -b ${bam_files} \
      -o ${DIR}/4_Images/${sv[3]}_${sv[1]}.png \
      -c ${sv[0]} \
      -s ${sv[1]} \
      -e ${e1} \
      -t ${sv[3]} \
      -w 200
  e2=$(( ${sv[2]} + ${sv[4]} ))
  echo ${sv[2]}
  echo $e2
  samplot.py \
      -n ${samples} \
      -b ${bam_files} \
      -o ${DIR}/4_Images/${sv[3]}_${sv[1]}_origin.png \
      -c ${sv[0]} \
      -s ${sv[2]} \
      -e ${e2} \
      -t ${sv[3]} \
      -d 300

  echo Done ${sv[0]} ${sv[1]} ${sv[3]}
  IFS=$'\n'
done
