#!/bin/bash

# Create Agg CSV for Cellranger Agg

# Assign arguments to variables
OUT_DIR=$1
INPUT_FILE=$2
OUTPUT_FILE=$3

# Add header 
echo "sample_id,molecule_h5" > "${OUTPUT_FILE}"

while IFS= read -r line
do
    echo "${line},${OUT_DIR}/GEX/${line}/outs/molecule_info.h5" >> "${OUTPUT_FILE}"
done < "${INPUT_FILE}"
