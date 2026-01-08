.libPaths('/hpc/group/dhvi/WieheLab/scripts/nextflow_pipelines/10x-nextflow-pipeline/bin/R_Libs')
library(ggplot2)
library(gtools)
library(dplyr)
library(stringr)
library(data.table)
library(qdap)
library(seqinr)
library(plotly)
library(Cairo)

args <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]
paths_file <- args[2]
feature_reference <- args[3]
vdj_path <- args[4]
Exp <- args[5]
PTID <- args[6]
beam_control <- args[7]

mat_paths <- readLines(paths_file)

# setwd('/datacommons/dhvi/10X_Genomics_data/COVID_Pediatric_BEAM_06122024/S0016B1/D1287/CR7.2/')
# 
# outdir <- '/datacommons/dhvi/10X_Genomics_data/COVID_Pediatric_BEAM_06122024/S0016B1/D1287/CR7.2/BEAM/'
# paths_file <- 'BEAM/separated_paths.txt'
# feature_reference <- 'BEAM/cellranger/feature_reference.csv'
# vdj_path <- 'VDJ/'
# Exp <- 'COVIDPediatricBEAM'
# PTID <- 'S0016B1day1287'
# beam_control <- 'HSA'
# 
# mat_paths <- readLines(paths_file)


# Function to extract lane number and read the csv file
read_and_process_csv <- function(path) {
  lane_number <- str_extract(path, "\\d+(?=\\.mat2csv\\.out\\.csv)")
  print(lane_number)
  data <- fread(path)
  colnames(data) <- gsub('-1', paste0('-', lane_number), colnames(data))
  colnames(data)[1] <- "Antigen"
  data <- melt(data, id.vars = 'Antigen', variable.name = 'barcode', value.name = 'umis')
  data <- dcast(data, formula = barcode ~ Antigen, value.var = 'umis')
  return(data)
}

# Read and process each file
mat2csv_out_list <- lapply(mat_paths, read_and_process_csv)

# Combine all data
#combined_data <- rbindlist(mat2csv_out_list)
combined_data <- smartbind(list = mat2csv_out_list)

# Read VDJ data
vdj <- fread(paste0(vdj_path, "/analysis/clonal_analysis/10x_merged_clones.fxnl_1to1.csv"))
colnames(vdj)[1] <- 'barcode'

# Control parameter
control <- beam_control
antigen_names <- read.csv(feature_reference)
#antigen_names[name == 'No_Antibody', name := 'No_Antigen']

all <- merge(vdj, combined_data, by = 'barcode', all.x = TRUE) 
colnames(all) <- ifelse(grepl("^[0-9]", colnames(all)), paste0("X", colnames(all)), colnames(all))


#S is the UMI count for the target antigen, and N is the UMI count for the negative control antigen.

SignalPRIOR <- 1
NoisePRIOR <- 3

#antigens <- grep('BEAM[01-16]', colnames(all), value = TRUE)
#antigens <- setdiff(antigens, control)
antigens <- colnames(all)[(ncol(vdj)+1):(ncol(vdj)+nrow(antigen_names))]
antigens <- antigens[!antigens == control]
antigens2 <- ifelse(grepl("^[0-9]", antigens), paste0("X", antigens), antigens)

calculate_beam_score <- function(antigen_umis, control_umis, SignalPRIOR = 1, NoisePRIOR = 3){
  score <- (1 - pbeta(0.925, antigen_umis + SignalPRIOR, control_umis + NoisePRIOR)) * 100
  return(score)
}

scores <- data.frame(cbind(barcode = all$barcode, apply(all[, ..antigens], 2, function(x) calculate_beam_score(antigen_umis = x, control_umis = all[[control]])))) # Change to control column
scores <- scores %>% mutate_at(antigens2, as.numeric)

all_with_scores <- merge(all, scores, by = 'barcode', suffixes = c('', '_BEAM_Score'), all = T) 

order <- select(all_with_scores, paste0(antigens, '_BEAM_Score')) %>% apply(., 1, sum) %>% order(., decreasing = T)

all_with_scores <- all_with_scores[order,]

colnames(all_with_scores) <- mgsub(pattern = antigen_names$id, replacement = paste0(antigen_names$name, '_UMIs'), text.var = colnames(all_with_scores)) 
colnames(all_with_scores) <- gsub('_UMIs_BEAM_Score', '_BEAM_Score', colnames(all_with_scores))
all_with_scores <- all_with_scores %>% mutate_at(vars(ends_with("BEAM_Score")), ~round(., 3))

#write.csv(all_with_scores, file = paste0(Exp, '_', PTID, '_merged_clonal_analysis.fxnl_1to1.BEAM_Scores.csv'), row.names = F, quote = F)

# Add sequences
#heavy_seqs <- read.fasta(paste0(vdj_path, "/analysis/unfiltered/Heavy/heavies.VDJ.fasta"), as.string = TRUE)
#kappa_seqs <- read.fasta(paste0(vdj_path, "/analysis/unfiltered/Kappa/kappas.VJ.fasta"), as.string = TRUE)
#lambda_seqs <- read.fasta(paste0(vdj_path, "/analysis/unfiltered/Lambda/lambdas.VJ.fasta"), as.string = TRUE)

#heavies_df <- data.frame(ReadID_Heavy = names(heavy_seqs), Heavy_VDJ_Seq = toupper(as.character(heavy_seqs)))
#kappas_df <- data.frame(ReadID_Light = names(kappa_seqs), Light_VJ_Seq = toupper(as.character(kappa_seqs)))
#lambdas_df <- data.frame(ReadID_Light = names(lambda_seqs), Light_VJ_Seq = toupper(as.character(lambda_seqs)))

#lights_df <- smartbind(kappas_df, lambdas_df)

#all_with_scores_and_seqs <- merge(all_with_scores, heavies_df, by = 'ReadID_Heavy', all.x = TRUE)
#all_with_scores_and_seqs <- merge(all_with_scores_and_seqs, lights_df, by = 'ReadID_Light', all.x = TRUE)

converted_barcodes <- read.csv('/datacommons/dhvi/mb488/scripts/10x/converted_10x_barcodes.csv', header = TRUE)

all_with_scores_and_seqs <- all_with_scores

all_with_scores_and_seqs$ID <- gsub('-[1-8]', '', all_with_scores_and_seqs$barcode)
all_with_scores_and_seqs$Gem_Well <- gsub('[AGCT]{16}-', '', all_with_scores_and_seqs$barcode)
all_with_scores_and_seqs$Light_chain_convert <- case_when(
  all_with_scores_and_seqs$Chain == 'Kappa' ~ 'K',
  all_with_scores_and_seqs$Chain == 'Lambda' ~ 'L'
)
all_with_scores_and_seqs$Heavy_contig <- gsub('[ACGT]{16}-[1-8]_contig_', '', all_with_scores_and_seqs$ReadID_Heavy)
all_with_scores_and_seqs$Light_contig <- gsub('[ACGT]{16}-[1-8]_contig_', '', all_with_scores_and_seqs$ReadID_Light)

all_with_scores_and_seqs <- merge(all_with_scores_and_seqs, converted_barcodes, by.x = 'ID', by.y = 'X10X_BC')

all_with_scores_and_seqs$ConvertedID <- paste(Exp, PTID, paste('G', all_with_scores_and_seqs$Gem_Well, sep = ''), all_with_scores_and_seqs$Converted_BC, sep = '_')
all_with_scores_and_seqs$Converted_ReadID_Heavy <- paste(Exp, PTID, paste('G', all_with_scores_and_seqs$Gem_Well, sep = ''), all_with_scores_and_seqs$Converted_BC, paste('H', all_with_scores_and_seqs$Heavy_contig, sep = ''), sep = '_')
all_with_scores_and_seqs$Converted_ReadID_Light <- paste(Exp, PTID, paste('G', all_with_scores_and_seqs$Gem_Well, sep = ''), all_with_scores_and_seqs$Converted_BC, paste(all_with_scores_and_seqs$Light_chain_convert, all_with_scores_and_seqs$Light_contig, sep = ''), sep = '_')

all_with_scores_and_seqs <- select(all_with_scores_and_seqs, ConvertedID, barcode, Fxnl_1to1, Converted_ReadID_Heavy, ReadID_Heavy, CloneID_Heavy, X.Members_Heavy, VGene_Heavy, 
                                   DGene, JGene_Heavy, CDR3Length_Heavy, Isotype, MuFreq_Heavy, Converted_ReadID_Light, ReadID_Light, CloneID_Light, X.Members_Light, Chain,
                                   VGene_Light, JGene_Light, CDR3Length_Light, MuFreq_Light, Heavy_Fxnl, Light_Fxnl, N_Heavy, N_Kappa, N_Lambda, N_Light, umis_Heavy, umis_Light, 
                                   CDR3_Heavy, CDR3_Light, CDR3_Heavy_AA, CDR3_Light_AA, Heavy_VDJ_Seq, Light_VJ_Seq, contains('UMIs'), contains('BEAM_Score'))

order <- select(all_with_scores_and_seqs, contains('_BEAM_Score')) %>% apply(., 1, sum) %>% order(., decreasing = T)

all_with_scores_and_seqs <- all_with_scores_and_seqs[order, ]
if ('DGene' %in% colnames(all_with_scores_and_seqs)){
  colnames(all_with_scores_and_seqs)[which(colnames(all_with_scores_and_seqs)=='DGene')] <- 'DGene_Heavy'
}

#colnames(all_with_scores_and_seqs)[colnames(all_with_scores_and_seqs) %in% antigen_names$name] <- paste(colnames(all_with_scores_and_seqs)[colnames(all_with_scores_and_seqs) %in% antigen_names$name], '_UMIs', sep = '')

write.csv(all_with_scores_and_seqs, file = paste0(outdir, '/', Exp, '_', PTID, '_merged_clonal_analysis.fxnl_1to1.BEAM_Scores.UMIs.Seqs.csv'), row.names = FALSE, quote = FALSE)

# Get long format for ggplot input
cols <- colnames(all_with_scores)[grepl('_BEAM_Score', colnames(all_with_scores))]

heatmap_sub <- select(all_with_scores, barcode, all_of(cols)) %>% data.frame()
heatmap_sub <- heatmap_sub[!is.na(heatmap_sub[,2]),]

heatmap_long <- reshape2::melt(heatmap_sub, id.vars=c("barcode"))
colnames(heatmap_long) <- c('Barcode', 'Antigen', 'BEAM_Score')

heatmap_long$Antigen <- gsub('_BEAM_Score', '', heatmap_long$Antigen)
antigen_levels <- ifelse(grepl("^[0-9]",antigen_names$name), paste0("X", antigen_names$name), antigen_names$name)
antigen_levels <- antigen_levels[-which(antigen_levels == beam_control)]

heatmap_long$Barcode <- factor(heatmap_long$Barcode, levels = all_with_scores$barcode)
heatmap_long$Antigen <- factor(heatmap_long$Antigen, levels = antigen_levels)

width <- nrow(heatmap_sub)/1000 + 6
height <- length(antigens) + 1

CairoPDF(paste0(outdir, '/', Exp, '_', PTID, ".heatmap.pdf"), height = height, width = width)
ggplot(heatmap_long, aes(x = Barcode, y = Antigen, fill = BEAM_Score)) + geom_tile() + 
  scale_fill_gradient(low = 'grey90', high = 'red', limits = c(0,100)) + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.title = element_text(hjust = .5)) + 
  xlab('Barcode') + 
  ylab('Antigen') + labs(fill = 'BEAM Score') + scale_y_discrete(expand = c(0,0)) + scale_x_discrete(expand = c(0,0))
dev.off()

# Get numeric cols
numeric_cols <- colnames(heatmap_sub)[-1]

# Apply sum function to numeric columns
heatmap_sub$rowsum <- apply(heatmap_sub[, numeric_cols], 1, sum)

# Subset based on condition that row sums is greater than 10 or > nrow(feature_reference)
nonzero_barcodes <- heatmap_sub[heatmap_sub$rowsum > nrow(antigen_names), 'barcode']

# Check if nonzero barcodes is not empty 
if (length(nonzero_barcodes) > 0 ) {
  heatmap_long_sub <- subset(heatmap_long, Barcode %in% nonzero_barcodes)

  heatmap_out_trimmed_name <- paste0(Exp, '_', PTID, ".heatmap.trimmed.pdf")
  
  width <- length(nonzero_barcodes)/1000 + 10
  height <- length(antigens) + 1

  CairoPDF(paste0(outdir, '/', heatmap_out_trimmed_name), height = height, width = width)
  h <- ggplot(heatmap_long_sub, aes(x = Barcode, y = Antigen, fill = BEAM_Score)) + geom_tile() + 
    scale_fill_gradient(low = 'grey90', high = 'red', limits = c(0,100)) + 
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.title = element_text(hjust = .5)) + 
    xlab('Barcode') + 
    ylab('Antigen') + labs(fill = 'BEAM Score') + scale_y_discrete(expand = c(0,0)) + scale_x_discrete(expand = c(0,0)) 
  print(h)
  dev.off()
} else {
  message("No nonzero barcodes found. Skipping trimmed heatmap generation.")
}

# heatmap_long$BEAM_Score <- round(heatmap_long$BEAM_Score, digits = 3)
# 
# plotly_out_trimmed_name <- paste0(Exp, '_', PTID, ".heatmap.long.html")
# 
# options(bitmapType="cairo")
# #pdf(plotly_out_trimmed_name, height = 8, width = 28)
# p <- ggplot(heatmap_long, aes(x = Barcode, y = Antigen, fill = BEAM_Score)) + geom_tile() + 
#   scale_fill_gradient(low = 'grey90', high = 'red', limits=c(0,100)) + 
#   theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.title = element_text(hjust = .5)) + 
#   xlab('Barcode') + 
#   ylab('Antigen') + labs(fill = 'BEAM Score') + scale_y_discrete(expand = c(0,0)) + scale_x_discrete(expand = c(0,0))
# 
# l <- plotly::ggplotly(p)
# 
# htmlwidgets::saveWidget(l, file = plotly_out_trimmed_name)
