# 10x Nextflow Pipeline

## Overview

This Nextflow pipeline is designed for processing 10X Genomics data, including Gene Expression (GEX), VDJ, and BEAM. There are three main workflows: GEX, VDJ, and BEAM. GEX uses Cellranger's mkfastq,count and aggr functions. VDJ uses Cellranger's mkfastq and VDJ functions, an in-house aggregation script, and in-house Cloanalyst scripts for analysis. BEAM uses Cellranger's mkfastq, count, and mat2csv functions, and an in-house script to calculate BEAM scores. The BEAM workflow incorporates both VDJ and BEAM, so there is no need to run these workflows separately when running BEAM. 

## Getting started

- Nextflow is already installed on the cluster.
- Make sure that you have either the Corretto or Temurin (LTS) distributions of Java installed. Other versions might work initially but can lead to difficult bugs over time.
    - Install Java through SKDMAN! (Software Development Kit Manager)
- Ensure that you are running the Nextflow job script from the parent directory that you would like your output to be stored in. 
- You must run the job script wrapper (Ex. sbatch nextflow_gex.job) Be sure to modify the arguments in the script. (Will eventually update this to pass arguments to the job script itself)
- There is a template job script for each workflow (nextflow_gex.job | nextflow_vdj.job | nextflow_beam.job)
- /cwork/[USER_ID]/work will automatically be set as the Nextflow working directory to store each process's log files, intermediate files, and outputs.
- The Nextflow process log will be automatically created in the directory that you run Nextflow from.

## Parameters
* exp: Experiment name (Only required for BEAM)
* ptid: Patient ID (Only required for BEAM. If you need to run GEX separately for different lanes when using different libraries, specifiy the lane at the end of the PTID, eg. "N_2976_507_lane1")
* sample_id: Sample name in the CSV (ex: If sample column in CSV has COVID_BEAM_1, COVID_BEAM_2 etc, the sample_id will be COVID_BEAM).
* species: 'human' | 'rhesus' | 'mouse'
* gex_raw: Path to the sequencing run directory for GEX.
* vdj_raw: Path to the sequencing run directory for VDJ.
* beam_raw: Path to the sequencing run directory for BEAM.
* gex_csv: Path to the required csv for the GEX workflow.
* vdj_csv: Path to the required csv for the VDJ workflow.
* beam_csv: Path to the required csv for the BEAM workflow.
* workflow: 'gex' | 'vdj' | 'beam'
* outdir: Output directory (Parent directory to house the output folders. The pipeline will automatically generate folders called "VDJ", "GEX", and "BEAM", so ensure that the output directory is one level above these.)
* cellranger: Version of cellranger to use
* feature_ref: Feature reference CSV. **(This is only required for the BEAM job script)**
* clonal: 'c' | 'p' | 'b'  Cloanalyst, partis, or both

Note: The transcriptome and vdj references for human, rhesus, and mouse are already pre-set and will be used accordingly depending on the species argument specified. You can change this at any time by adding these paramaters to the job script:
* --human_transcriptome [path]
* --rhesus_transcriptome [path]
* --mouse_transcriptome [path]
* --human_vdj_ref [path]
* --rhesus_vdj_ref [path]
* --mouse_vdj_ref [path]

## Usage
```
sbatch nextflow_gex.job
sbatch nextflow_vdj.job
sbatch nextflow_beam.job
```

## Resources
* Nextflow documentation: https://www.nextflow.io/docs/latest/index.html
* Cellranger Documentation: https://www.nextflow.io/docs/latest/index.html
