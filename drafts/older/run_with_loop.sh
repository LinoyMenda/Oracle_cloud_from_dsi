#!/bin/bash

# Define the path to the accession file
SRR_ACC_FILE_PATH="./top10_SRR.txt"

# Read the accession numbers from the file into an array
mapfile -t srr_acc_ids < "$SRR_ACC_FILE_PATH"

# Loop over the accession numbers and download each one using fasterq-dump
for srr_id in "${srr_acc_ids[@]}"; do
    # Do something with each accession number, e.g., download using fasterq-dump
    echo "Processing accession: $srr_id"
    python src/CDR3_extract_algo.py "${srr_id}" -d true

    #if [[ -f "./algo_results/CDR3_results_${srr_id}.txt" ]]; then
        #oci os object put --bucket-name algo_linoy --file "./algo_results/CDR3_results_${srr_id}.txt"
   #else
        #echo "nothing to push"
   #fi
done

