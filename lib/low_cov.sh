#!/bin/bash

while getopts ":s:o:" flag; do
  case "${flag}" in
    s ) SAMPLE="${OPTARG}" ;;
    o ) OUTDIR="${OPTARG}" ;;
  esac
done
