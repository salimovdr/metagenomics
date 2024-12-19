#!/bin/bash

INPUT_FILE="data/scaffolds_TB2_SN3.fasta"
OUTPUT_DIR="maxbin_results"
BASENAME="TB2_SN3"

mkdir -p $OUTPUT_DIR

grep ">" "$INPUT_FILE" | cut -f 2 -d ">" > "${BASENAME}_names.txt"
grep ">" "$INPUT_FILE" | cut -f 6 -d "_" > "${BASENAME}_cov.txt"
paste "${BASENAME}_names.txt" "${BASENAME}_cov.txt" > "${BASENAME}_abund.txt"

run_MaxBin.pl -contig "$INPUT_FILE" -abund "${BASENAME}_abund.txt" -out "${OUTPUT_DIR}/${BASENAME}" -thread 15

rm "${BASENAME}_names.txt" "${BASENAME}_cov.txt" "${BASENAME}_abund.txt"

echo "MaxBin2 завершён. Результаты находятся в папке ${OUTPUT_DIR}/"
