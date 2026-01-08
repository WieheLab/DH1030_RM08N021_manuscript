#!/bin/bash

conda activate /datacommons/dhvi/srb108/partis
module load R/4.1.1-rhel8

# Get folder with CellRanger VDJ output
LANE_FOLDER=$1

# Get the current user
USERNAME=$(id -un)

Rscript /hpc/group/dhvi/WieheLab/scripts/nextflow_pipelines/10x-nextflow-pipeline/bin/get_vdj_all_lane_files.R ${LANE_FOLDERS}

# Partis working directory
mkdir -p ./partis_workdir/
PARTIS_WORKDIR=$(readlink -f ./partis_workdir)


# Check if Partis working directory exists
# Check if the directory exists
if [ ! -d "$PARTIS_WORKDIR" ]; then
  # If it doesn't exist, create it
  mkdir -p "$PARTIS_WORKDIR"
fi

# Export the Partis working directory environment variable
export PARTIS_WORKDIR
echo "The Partis working directory is set to $PARTIS_WORKDIR"

mkdir partis_analysis
cp all_contig_annotations.all_lanes.csv ./partis_analysis  
cp all_contig.all_lanes.fasta ./partis_analysis
cd partis_analysis
mkdir unfiltered cellranger_output filtered intermediate_files clonal_analysis

# Move cellranger files and split cellranger output by chain for input into Partis
mv all_contig* cellranger_output
cd cellranger_output

Rscript /datacommons/dhvi/mb488/scripts/10x/R/all_contig_annots_umi_filter.R -a all_contig_annotations.all_lanes.csv -o all_contig_annotations.all_lanes.umi_filtered.csv
python /datacommons/dhvi/mb488/scripts/10x/python/split_fasta_by_chain_filter_cell_high_conf_CR6.py all_contig_annotations.all_lanes.umi_filtered.csv all_contig.all_lanes.fasta 


# Create folders for Partis and run each chain
cd ../unfiltered
mkdir Heavy Kappa Lambda
mv ../cellranger_output/heavies.fasta Heavy
mv ../cellranger_output/kappas.fasta Kappa 
mv ../cellranger_output/lambdas.fasta Lambda

cd Heavy 
# First check if the output directory exists and remove it if exists. Output directory needs to be deleted on each re-run.
DIR_NAME="_partis_output_igh"

# Check if the directory exists
if [ -d "$DIR_NAME" ]; then
  # If it exists, delete it
  rm -rf "$DIR_NAME"
  echo "Directory $DIR_NAME has been deleted."
else
  echo "Directory $DIR_NAME does not exist."
fi

# Annotate with Partis
job_id=$(sbatch /hpc/group/dhvi/WieheLab/scripts/partis/job_scripts/partis_flexible_annotate_passParamDir.job \
-p heavies.fasta \
-c igh \
-g /datacommons/dhvi/partis_germlines/KimDb_Rhesus_wo_tryp_pos \
-w $PARTIS_WORKDIR \
-m /datacommons/dhvi/Projects/DHVI/Wiehe/Partis_Hash_Create_2024/partis_rhesus10X/heavy_kimdb/_partis_output_igh \
-t 10 \
-v | awk '{print $4}')

echo "Submitted partis job for heavy chain annotation with ID: $job_id"

# Wait for the job to finish
while true; do
    # Get the job state using scontrol
    job_state=$(scontrol show job $job_id | grep JobState | awk -F= '{print $2}' | awk '{print $1}')
    
    # Check if the job is completed
    if [[ "$job_state" == "COMPLETED" ]]; then
        echo "Job $job_id completed successfully."
        break
    elif [[ "$job_state" == "FAILED" || "$job_state" == "CANCELLED" ]]; then
        echo "Job $job_id failed or was cancelled."
        exit 1
    else
        # Wait for a bit before checking again
        sleep 10
    fi
done


# Convert to SMUA, Recombination Summaries Files
python3 /hpc/group/dhvi/WieheLab/scripts/partis/python_scripts/partis_to_RS_IMGT_hashed_wErr.py -p heavies.partis.annotate.yaml -c heavy -s other -a /datacommons/dhvi/partis_germlines/IMGT_V_region_annotation_HASHES/kimdb_auto_update.manual4.yaml -r 1 -n 0
# VDJ trimming
Rscript /hpc/group/dhvi/WieheLab/scripts/R/scripts/vdjtrimmer.R -s VH_SimpleMarkedUAs.fasta -o heavies.VDJ.fasta

HEAVY_FASTA=$(readlink -f heavies.fasta)

cd ../Kappa 
# First check if the output directory exists and remove it if exists. Output directory needs to be deleted on each re-run.
DIR_NAME="_partis_output_igk"

# Check if the directory exists
if [ -d "$DIR_NAME" ]; then
  # If it exists, delete it
  rm -rf "$DIR_NAME"
  echo "Directory $DIR_NAME has been deleted."
else
  echo "Directory $DIR_NAME does not exist."
fi

# Annotate with partis
job_id=$(sbatch /hpc/group/dhvi/WieheLab/scripts/partis/job_scripts/partis_flexible_annotate_passParamDir.job \
-p kappas.fasta \
-c igk \
-g /datacommons/dhvi/Projects/DHVI/Wiehe/Partis_Hash_Create_2024/partis_rhesus10X/BW25_preprint/kappa/_partis_output_igk/hmm/germline-sets \
-w $PARTIS_WORKDIR \
-m /datacommons/dhvi/Projects/DHVI/Wiehe/Partis_Hash_Create_2024/partis_rhesus10X/BW25_preprint/kappa/_partis_output_igk \
-t 10 \
-v | awk '{print $4}')

echo "Submitted partis job for kappa annotation with ID: $job_id"

# Wait for the job to finish
while true; do
    # Get the job state using scontrol
    job_state=$(scontrol show job $job_id | grep JobState | awk -F= '{print $2}' | awk '{print $1}')
    
    # Check if the job is completed
    if [[ "$job_state" == "COMPLETED" ]]; then
        echo "Job $job_id completed successfully."
        break
    elif [[ "$job_state" == "FAILED" || "$job_state" == "CANCELLED" ]]; then
        echo "Job $job_id failed or was cancelled."
        exit 1
    else
        # Wait for a bit before checking again
        sleep 10
    fi
done

# Convert to SMUA, Recombination Summaries Files
python3 /hpc/group/dhvi/WieheLab/scripts/partis/python_scripts/partis_to_RS_IMGT_hashed_wErr.py -p kappas.partis.annotate.yaml -c kappa -s rhesus -r 1 -n 0 

# VDJ trimming
Rscript /hpc/group/dhvi/WieheLab/scripts/R/scripts/vdjtrimmer.R -s VK_SimpleMarkedUAs.fasta -o kappas.VJ.fasta

KAPPA_FASTA=$(readlink -f kappas.fasta)
#PARTIS_KAPPA_OUT=$(readlink -f _partis_output_igk)


cd ../Lambda 
# First check if the output directory exists and remove it if exists. Output directory needs to be deleted on each re-run.
DIR_NAME="_partis_output_igl"

# Check if the directory exists
if [ -d "$DIR_NAME" ]; then
  # If it exists, delete it
  rm -rf "$DIR_NAME"
  echo "Directory $DIR_NAME has been deleted."
else
  echo "Directory $DIR_NAME does not exist."
fi

# Annotate with partis
job_id=$(sbatch /hpc/group/dhvi/WieheLab/scripts/partis/job_scripts/partis_flexible_annotate_passParamDir.job \
-p lambdas.fasta \
-c igl \
-g /datacommons/dhvi/Projects/DHVI/Wiehe/Partis_Hash_Create_2024/partis_rhesus10X/BW25_preprint/lambda/_partis_output_igl/hmm/germline-sets \
-w $PARTIS_WORKDIR \
-m /datacommons/dhvi/Projects/DHVI/Wiehe/Partis_Hash_Create_2024/partis_rhesus10X/BW25_preprint/lambda/_partis_output_igl \
-t 10 \
-v | awk '{print $4}')

echo "Submitted partis job for lambda annotation with ID: $job_id"

# Wait for the job to finish
while true; do
    # Get the job state using scontrol
    job_state=$(scontrol show job $job_id | grep JobState | awk -F= '{print $2}' | awk '{print $1}')
    
    # Check if the job is completed
    if [[ "$job_state" == "COMPLETED" ]]; then
        echo "Job $job_id completed successfully."
        break
    elif [[ "$job_state" == "FAILED" || "$job_state" == "CANCELLED" ]]; then
        echo "Job $job_id failed or was cancelled."
        exit 1
    else
        # Wait for a bit before checking again
        sleep 10
    fi
done



# Convert to SMUA, Recombination Summaries Files
python3 /hpc/group/dhvi/WieheLab/scripts/partis/python_scripts/partis_to_RS_IMGT_hashed_wErr.py -p lambdas.partis.annotate.yaml -c lambda -s rhesus -r 1 -n 0 

# VDJ trimming
Rscript /hpc/group/dhvi/WieheLab/scripts/R/scripts/vdjtrimmer.R -s VL_SimpleMarkedUAs.fasta -o lambdas.VJ.fasta

LAMBDA_FASTA=$(readlink -f lambdas.fasta)
#PARTIS_LAMBDA_OUT=$(readlink -f _partis_output_igl)


# Get chain count for each chain
cd ../../intermediate_files

Rscript /datacommons/dhvi/mb488/scripts/10x/R/get_chain_counts.R -h ../unfiltered/Heavy/VH_RecombinationSummaries.RF.txt -k ../unfiltered/Kappa/VK_RecombinationSummaries.RF.txt -l ../unfiltered/Lambda/VL_RecombinationSummaries.RF.txt

# Get Heavy, Kappa, Lambda IDs
cut -f 1 ../unfiltered/Heavy/VH_RS.RF.fxnl.prod.txt > heavies.fxnl.ids.txt
cut -f 1 ../unfiltered/Kappa/VK_RS.RF.fxnl.prod.txt > kappas.fxnl.ids.txt
cut -f 1 ../unfiltered/Lambda/VL_RS.RF.fxnl.prod.txt > lambdas.fxnl.ids.txt

cd ../clonal_analysis

mkdir Heavy Kappa Lambda 

cd Heavy 
cp ../../unfiltered/Heavy/VH_RecombinationSummaries.RF.txt . 

# Partition with Partis
job_id=$(sbatch /hpc/group/dhvi/WieheLab/scripts/partis/job_scripts/partis_flexible_partition_passParamDir.job \
-p $HEAVY_FASTA \
-c igh \
-g /datacommons/dhvi/partis_germlines/KimDb_Rhesus_wo_tryp_pos \
-w $PARTIS_WORKDIR \
-t 10 \
-m /datacommons/dhvi/Projects/DHVI/Wiehe/Partis_Hash_Create_2024/partis_rhesus10X/heavy_kimdb/_partis_output_igh | awk '{print $4}')

echo "Submitted partis job for heavy chain partitioning with ID: $job_id"

# Wait for the job to finish
while true; do
    # Get the job state using scontrol
    job_state=$(scontrol show job $job_id | grep JobState | awk -F= '{print $2}' | awk '{print $1}')
    
    # Check if the job is completed
    if [[ "$job_state" == "COMPLETED" ]]; then
        echo "Job $job_id completed successfully."
        break
    elif [[ "$job_state" == "FAILED" || "$job_state" == "CANCELLED" ]]; then
        echo "Job $job_id failed or was cancelled."
        exit 1
    else
        # Wait for a bit before checking again
        sleep 10
    fi
done

 

cp ../../unfiltered/Heavy/heavies.partis.partition.yaml .
# Parse clone assignments with Partis
ABS_RF=$(readlink -f VH_RecombinationSummaries.RF.txt)
ABS_YML=$(readlink -f heavies.partis.partition.yaml)
python3 /hpc/group/dhvi/WieheLab/scripts/partis/python_scripts/get_clones.py -y $ABS_YML -r $ABS_RF 

cd ../Kappa  
cp ../../unfiltered/Kappa/VK_RecombinationSummaries.RF.txt .

# Partition with Partis
job_id=$(sbatch /hpc/group/dhvi/WieheLab/scripts/partis/job_scripts/partis_flexible_partition_passParamDir.job \
-p $KAPPA_FASTA \
-c igk \
-g /datacommons/dhvi/Projects/DHVI/Wiehe/Partis_Hash_Create_2024/partis_rhesus10X/BW25_preprint/kappa/_partis_output_igk/hmm/germline-sets \
-w $PARTIS_WORKDIR \
-t 10 \
-m /datacommons/dhvi/Projects/DHVI/Wiehe/Partis_Hash_Create_2024/partis_rhesus10X/BW25_preprint/kappa/_partis_output_igk | awk '{print $4}')

echo "Submitted partis job for kappa partitioning with ID: $job_id"

# Wait for the job to finish
while true; do
    # Get the job state using scontrol
    job_state=$(scontrol show job $job_id | grep JobState | awk -F= '{print $2}' | awk '{print $1}')
    
    # Check if the job is completed
    if [[ "$job_state" == "COMPLETED" ]]; then
        echo "Job $job_id completed successfully."
        break
    elif [[ "$job_state" == "FAILED" || "$job_state" == "CANCELLED" ]]; then
        echo "Job $job_id failed or was cancelled."
        exit 1
    else
        # Wait for a bit before checking again
        sleep 10
    fi
done

cp ../../unfiltered/Kappa/kappas.partis.partition.yaml .
# Parse clone assignments with Partis
ABS_RF=$(readlink -f VK_RecombinationSummaries.RF.txt)
ABS_YML=$(readlink -f kappas.partis.partition.yaml)
python3 /hpc/group/dhvi/WieheLab/scripts/partis/python_scripts/get_clones.py -y $ABS_YML -r $ABS_RF 

cd ../Lambda 

cp ../../unfiltered/Lambda/VL_RecombinationSummaries.RF.txt .

# Partition with Partis
job_id=$(sbatch /hpc/group/dhvi/WieheLab/scripts/partis/job_scripts/partis_flexible_partition_passParamDir.job \
-p $LAMBDA_FASTA \
-c igl \
-g /datacommons/dhvi/Projects/DHVI/Wiehe/Partis_Hash_Create_2024/partis_rhesus10X/BW25_preprint/lambda/_partis_output_igl/hmm/germline-sets \
-w $PARTIS_WORKDIR \
-t 10 \
-m /datacommons/dhvi/Projects/DHVI/Wiehe/Partis_Hash_Create_2024/partis_rhesus10X/BW25_preprint/lambda/_partis_output_igl | awk '{print $4}')

echo "Submitted partis job for lambda partitioning with ID: $job_id"

# Wait for the job to finish
while true; do
    # Get the job state using scontrol
    job_state=$(scontrol show job $job_id | grep JobState | awk -F= '{print $2}' | awk '{print $1}')
    
    # Check if the job is completed
    if [[ "$job_state" == "COMPLETED" ]]; then
        echo "Job $job_id completed successfully."
        break
    elif [[ "$job_state" == "FAILED" || "$job_state" == "CANCELLED" ]]; then
        echo "Job $job_id failed or was cancelled."
        exit 1
    else
        # Wait for a bit before checking again
        sleep 10
    fi
done

cp ../../unfiltered/Lambda/lambdas.partis.partition.yaml .

# Parse clone assignments with Partis
ABS_RF=$(readlink -f VL_RecombinationSummaries.RF.txt)
ABS_YML=$(readlink -f lambdas.partis.partition.yaml)
python3 /hpc/group/dhvi/WieheLab/scripts/partis/python_scripts/get_clones.py -y $ABS_YML -r $ABS_RF 

cd ../
Rscript /hpc/group/dhvi/WieheLab/scripts/nextflow_pipelines/10x-nextflow-pipeline/bin/partis_merge_clone_info.R ../cellranger_output/all_contig_annotations.all_lanes.umi_filtered.csv

cd ../
Rscript /hpc/group/dhvi/WieheLab/scripts/nextflow_pipelines/srb108/10x-nextflow-pipeline/bin/partis_subset_rs_from_merged_clonal_analysis.R

cd ../ 