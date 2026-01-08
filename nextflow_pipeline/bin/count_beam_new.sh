#!/bin/bash

# Assign the arguments to variables
FASTQS_DIR=$1
SAMPLE_ID=$2
TRANSCRIPTOME=$3
CELLRANGER=$4
FEATURE_REF=$5

# Get the name of the folder from mkfastq_outs that has the fastqs
FASTQ_FOLDER=$(ls -d ${FASTQS_DIR}/*/ | grep -vE "${FASTQS_DIR}/(\\.nextflow|Reports|Stats)/$" | awk -F '/' '{print $(NF-1)}')
echo "FASTQ_FOLDER: ${FASTQ_FOLDER}"

touch library.csv
echo "fastqs,sample,library_type" > library.csv
echo "${FASTQS_DIR}/${FASTQ_FOLDER},${SAMPLE_ID},Antibody Capture" >> library.csv

# Run cellranger count
/datacommons/dhvi/10X_Genomics_data/cellranger/cellranger-${CELLRANGER}/cellranger count --id=${SAMPLE_ID} --transcriptome=${TRANSCRIPTOME} --libraries=library.csv --feature-ref=${FEATURE_REF} --create-bam false