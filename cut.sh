#!/bin/bash

output_dir="data/fastqs/trimmed_quality"

mkdir -p ${output_dir}

for file in data/fastqs/*.fastq.gz; do
    echo "Processing ${file}..."

    base=$(basename ${file})

    trimmomatic SE -phred33 ${file} ${output_dir}/trimmed_${base} SLIDINGWINDOW:4:20 MINLEN:50

    echo "${file} has been trimmed and saved to ${output_dir}/trimmed_${base}"
done

echo "All files have been processed. Trimmed files are in the ${output_dir} directory."
