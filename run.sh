#!/bin/bash
export PATH=$PATH:/root/bin
FILE_PATH="/app/algo_results"

python CDR3_extract_algo.py "${FILE_ID}"

if [ -f "${FILE_PATH}/CDR3_results_${FILE_ID}.txt" ]; then
    echo "Try to upload file: ${FILE_PATH}/CDR3_results_${FILE_ID}.txt ..."
    oci os object put --bucket-name bucket_Linoy --file "${FILE_PATH}/CDR3_results_${FILE_ID}.txt" --region us-ashburn-1
    echo "Check if exist in bucket ..."
else
    echo "nothing to push"
fi
