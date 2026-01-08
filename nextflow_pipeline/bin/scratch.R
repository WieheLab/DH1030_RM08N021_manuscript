#!/usr/bin/env Rscript
.libPaths('/datacommons/dhvi/mb488/R_Libs/R_v4.1.1/')

library(seqinr)

args <- commandArgs(trailingOnly = TRUE)
lane_folders <- args[1]

# Read the file containing directory paths and split by comma
paths_file <- readLines(lane_folders)  # Ensure this is the correct path to your file
directories <- unlist(strsplit(paths_file, split = ",", fixed = TRUE))

# Initialize a vector to store all file paths that match the pattern
all_files <- character()

# Loop through each directory to find files matching the pattern
for (dir in directories) {
  # Use list.files to find all files matching the pattern recursively
  files_in_dir <- list.files(path = dir, pattern = 'testing.txt$', recursive = TRUE, full.names = TRUE)
  # Combine the found files into the all_files vector
  all_files <- c(all_files, files_in_dir)
}

# Now all_files contains paths to all files matching 'all_contig.fasta$' across the specified directories
print('Found these files:')
print(all_files)