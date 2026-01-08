#!/bin/bash

module load R/4.1.1-rhel8

LANE_FOLDERS=$1
CLOANALYST_PATH=/datacommons/dhvi/cloanalyst-linux-mpi/CloanalystPackage

Rscript /hpc/group/dhvi/WieheLab/scripts/nextflow_pipelines/10x-nextflow-pipeline/bin/get_vdj_all_lane_files.R ${LANE_FOLDERS}

mkdir analysis
cp all_contig_annotations.all_lanes.csv ./analysis 
cp all_contig.all_lanes.fasta ./analysis
cd analysis

# Create directories
mkdir unfiltered cellranger_output filtered intermediate_files clonal_analysis

# Move cellranger files and split cellranger output by chain for input into Cloanalyst
mv all_contig* cellranger_output
cd cellranger_output

Rscript /datacommons/dhvi/mb488/scripts/10x/R/all_contig_annots_umi_filter.R -a all_contig_annotations.all_lanes.csv -o all_contig_annotations.all_lanes.umi_filtered.csv
python /datacommons/dhvi/mb488/scripts/10x/python/split_fasta_by_chain_filter_cell_high_conf_CR6.py all_contig_annotations.all_lanes.umi_filtered.csv all_contig.all_lanes.fasta 

# Create folders for Cloanalyst and run each chain
cd ../unfiltered
mkdir Heavy Kappa Lambda
mv ../cellranger_output/heavies.fasta Heavy
mv ../cellranger_output/kappas.fasta Kappa 
mv ../cellranger_output/lambdas.fasta Lambda

count_sequences() {
    grep -c "^>" "$1"
}

cd Heavy 
num_sequences=$(count_sequences heavies.fasta)
if [ "$num_sequences" -lt 10 ]; then
        num_threads=$num_sequences
    else
        num_threads=$SLURM_NTASKS_PER_NODE
    fi
echo "Using $num_threads threads for heavies.fasta"

#time mpirun -np $num_threads /opt/apps/rhel7/cloanalyst/cloanalyst_mpi-new parse --thread $num_threads -s "\"Mus musculus\"" -c heavy --excl 5 -g "\"AR20170307 FPC-F\"" heavies.fasta
mono ${CLOANALYST_PATH}/Cloanalyst.dll parse  --thread $num_threads -s "\"Mus musculus\"" -c heavy --excl 5 -g "\"AR20170307 FPC-F\"" heavies.fasta
/hpc/group/dhvi/WieheLab/scripts/Compiled/VDJ_trimmer/VDJ_trimmer -SMUA heavies.SimpleMarkedUAs.fasta -species mouse -chain heavy > heavies.VDJ.fasta

cd ../Kappa
num_sequences=$(count_sequences kappas.fasta)
if [ "$num_sequences" -lt 10 ]; then
        num_threads=$num_sequences
    else
        num_threads=$SLURM_NTASKS_PER_NODE
    fi
echo "Using $num_threads threads for kappas.fasta" 

#time mpirun -np $num_threads /opt/apps/rhel7/cloanalyst/cloanalyst_mpi-new parse --thread $num_threads -s "\"Mus musculus\"" -c kappa --excl 5 -g "\"AR20170307 FPC-F\"" kappas.fasta
mono ${CLOANALYST_PATH}/Cloanalyst.dll parse  --thread $num_threads -s "\"Mus musculus\"" -c kappa --excl 5 -g "\"AR20170307 FPC-F\"" kappas.fasta
/hpc/group/dhvi/WieheLab/scripts/Compiled/VDJ_trimmer/VDJ_trimmer -SMUA kappas.SimpleMarkedUAs.fasta -species mouse -chain kappa > kappas.VJ.fasta

cd ../Lambda 
num_sequences=$(count_sequences lambdas.fasta)
if [ "$num_sequences" -lt 10 ]; then
        num_threads=$num_sequences
    else
        num_threads=$SLURM_NTASKS_PER_NODE
    fi
echo "Using $num_threads threads for lambdas.fasta"

#time mpirun -np $num_threads /opt/apps/rhel7/cloanalyst/cloanalyst_mpi-new parse --thread $num_threads -s "\"Mus musculus\"" -c lambda --excl 5 -g "\"AR20170307 FPC-F\"" lambdas.fasta
mono ${CLOANALYST_PATH}/Cloanalyst.dll parse  --thread $num_threads -s "\"Mus musculus\"" -c lambda --excl 5 -g "\"AR20170307 FPC-F\"" lambdas.fasta
/hpc/group/dhvi/WieheLab/scripts/Compiled/VDJ_trimmer/VDJ_trimmer -SMUA lambdas.SimpleMarkedUAs.fasta -species	mouse -chain lambda > lambdas.VJ.fasta

# get chain count and fxnl calls for each chain
cd ../../intermediate_files

Rscript /datacommons/dhvi/mb488/scripts/10x/R/get_chain_counts.R -h ../unfiltered/Heavy/heavies.RecombinationSummaries.txt -k ../unfiltered/Kappa/kappas.RecombinationSummaries.txt -l ../unfiltered/Lambda/lambdas.RecombinationSummaries.txt

/hpc/group/dhvi/WieheLab/scripts/Compiled/RS_filter_functional/RS_filter_functional -c H -RS ../unfiltered/Heavy/heavies.RecombinationSummaries.txt -SMUA ../unfiltered/Heavy/heavies.SimpleMarkedUAs.fasta | cut -f1 > heavies.fxnl.ids.txt

/hpc/group/dhvi/WieheLab/scripts/Compiled/RS_filter_functional/RS_filter_functional -c K -RS ../unfiltered/Kappa/kappas.RecombinationSummaries.txt -SMUA ../unfiltered/Kappa/kappas.SimpleMarkedUAs.fasta | cut -f1 > kappas.fxnl.ids.txt

/hpc/group/dhvi/WieheLab/scripts/Compiled/RS_filter_functional/RS_filter_functional -c L -RS ../unfiltered/Lambda/lambdas.RecombinationSummaries.txt -SMUA ../unfiltered/Lambda/lambdas.SimpleMarkedUAs.fasta | cut -f1 > lambdas.fxnl.ids.txt

cd ../clonal_analysis

mkdir Heavy Kappa Lambda 

cd Heavy 
cp ../../unfiltered/Heavy/heavies.RecombinationSummaries.txt . 
mono $CLOANALYST_PATH/Cloanalyst.dll fast_partition heavies.RecombinationSummaries.txt

cd ../Kappa 
cp ../../unfiltered/Kappa/kappas.RecombinationSummaries.txt . 
mono $CLOANALYST_PATH/Cloanalyst.dll fast_partition kappas.RecombinationSummaries.txt

cd ../Lambda 
cp ../../unfiltered/Lambda/lambdas.RecombinationSummaries.txt . 
mono $CLOANALYST_PATH/Cloanalyst.dll fast_partition lambdas.RecombinationSummaries.txt

cd ../
Rscript /datacommons/dhvi/mb488/scripts/10x/R/merge_clone_info.R -I ../cellranger_output/all_contig_annotations.all_lanes.umi_filtered.csv

cd ../
Rscript /datacommons/dhvi/mb488/scripts/10x/R/subset_rs_from_merged_clonal_analysis.R