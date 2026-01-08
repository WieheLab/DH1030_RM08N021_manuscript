#!/bin/bash

# Assign the arguments to variables
FASTQS_DIR=$1
SAMPLE_ID=$2
TRANSCRIPTOME=$3
CELLRANGER=$4

# Get the name of the folder from mkfastq_outs that has the fastqs
FASTQ_FOLDER=$(ls -d ${FASTQS_DIR}/*/ | grep -vE "${FASTQS_DIR}/(\\.nextflow|Reports|Stats)/$" | awk -F '/' '{print $(NF-1)}')
echo "FASTQ_FOLDER: ${FASTQ_FOLDER}"

# Run cellranger count
/datacommons/dhvi/10X_Genomics_data/cellranger/cellranger-${CELLRANGER}/cellranger count --id=${SAMPLE_ID} --transcriptome=${TRANSCRIPTOME} --fastqs=${FASTQS_DIR}/${FASTQ_FOLDER} --sample=${SAMPLE_ID} --nosecondary --no-bam