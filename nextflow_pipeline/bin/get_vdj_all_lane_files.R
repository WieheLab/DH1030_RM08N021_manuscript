#.libPaths('/hpc/group/dhvi/WieheLab/scripts/nextflow_pipelines/10x-nextflow-pipeline/bin/R_Libs')

library(seqinr)
library(R.utils)  # For gunzip()

args <- commandArgs(trailingOnly = TRUE)
lane_folders <- args[1]

# Read the file containing directory paths
directories <- readLines(lane_folders)
directories <- directories[directories != "" & !grepl("^work$", directories)]

# Print the resulting list
print("Found these directories:")
print(directories)

# Initialize a vector to store all file paths that match the pattern
all_files <- character()

# Loop through each directory to find files matching the pattern
for (dir in directories) {
  # Use list.files to find all files matching the pattern recursively
  files_in_dir <- list.files(path = dir, pattern = 'all_contig.fasta$', recursive = TRUE, full.names = TRUE)
  # Combine the found files into the all_files vector
  all_files <- c(all_files, files_in_dir)
}

# Now all_files contains paths to all files matching 'all_contig.fasta$' across the specified directories
print('Found these files:')
print(all_files)

all_fasta <- data.frame()
all_annots <- data.frame()


for (i in all_files){
  print("printing file")
  print(i)

  # get lane number to substitute in fasta and annotations
  lane_num <- strsplit(i, split = '/')[[1]]
  lane_num <- lane_num[length(lane_num) - 2]  # Get second to last element
  
  lane_num <- unlist(strsplit(x = lane_num, split = '_'))
  lane_num <- lane_num[length(lane_num)]

  print("printing lane number")
  print(lane_num)
  
  # read in sequences and substitute lane number and create dataframe
  seqs <- read.fasta(i, forceDNAtolower = F, as.string = T)
  names(seqs) <- gsub('-1', paste0('-', lane_num), names(seqs))
  seqs_add <- data.frame(ID = names(seqs), seq = as.character(seqs))
  
  # read in annotations file and substute lane number in barcode and contig ID
  annots_file <- gsub('all_contig.fasta', 'all_contig_annotations.csv', i)
  annots <- read.csv(annots_file, header = T)
  annots$barcode <- gsub('-1', paste0('-', lane_num), annots$barcode)
  annots$contig_id <- gsub('-1', paste0('-', lane_num), annots$contig_id)
  
  # add lane to dataframe with all lanes and all sequences
  all_fasta <- rbind(all_fasta, seqs_add)
  all_annots <- rbind(all_annots, annots)

}

write.fasta(sequences = as.list(all_fasta$seq), names = all_fasta$ID, as.string = T, file.out = 'all_contig.all_lanes.fasta', nbchar = 100000)
write.csv(x = all_annots, file = 'all_contig_annotations.all_lanes.csv', quote = F, row.names = F)




