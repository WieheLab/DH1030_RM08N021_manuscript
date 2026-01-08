#!/bin/bash

# Assign the arguments as variables
FASTQS_DIR=$1
SAMPLE_ID=$2
REFERENCE=$3
CELLRANGER=$4
INNER_ENRICHMENT_PRIMERS_PATH=$5

# Get the name of the folder from mkfastq_outs that has the fastqs
FASTQ_FOLDER=$(ls -d ${FASTQS_DIR}/*/ | grep -vE "${FASTQS_DIR}/(\\.nextflow|Reports|Stats)/$" | awk -F '/' '{print $(NF-1)}')
echo "FASTQ_FOLDER: ${FASTQ_FOLDER}"

# Run cellranger vdj
/datacommons/dhvi/10X_Genomics_data/cellranger/cellranger-${CELLRANGER}/cellranger vdj --id=${SAMPLE_ID} --reference=${REFERENCE} --fastqs=${FASTQS_DIR}/${FASTQ_FOLDER} --sample=${SAMPLE_ID} --inner-enrichment-primers=${INNER_ENRICHMENT_PRIMERS_PATH} --chain IG
