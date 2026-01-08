#!/usr/bin/env Rscript
.libPaths(c('/hpc/group/dhvi/WieheLab/scripts/nextflow_pipelines/10x-nextflow-pipeline/bin/R_Libs'))

library(stringr)
library(gtools)

heavy_rs <- read.csv('../unfiltered/Heavy/VH_RecombinationSummaries.RF.txt', sep = "\t", header = T, stringsAsFactors = F)

heavy_rs$ID <- gsub('_contig_[0-9]{1,2}', '', heavy_rs$UID)
heavy_counts <- data.frame(table(heavy_rs$ID))
colnames(heavy_counts) <- c('ID', 'N_Heavy')

write.csv(heavy_counts, file="chain_counts.csv", quote = F, row.names = F)
