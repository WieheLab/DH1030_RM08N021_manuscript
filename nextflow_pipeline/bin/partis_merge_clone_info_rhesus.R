#.libPaths('/datacommons/dhvi/mb488/R_Libs/R_v4.0.0')

library(stringr)
library(dplyr)
library(optparse)
library(gtools)
library(Biostrings) 

args <- commandArgs(trailingOnly = TRUE)

contig_file <- args[1]

chain_counts <- read.csv("../intermediate_files/chain_counts.csv")
heavy_fxnl <- read.csv("../intermediate_files/heavies.fxnl.ids.txt", header = T, stringsAsFactors = F)
annots <- read.csv(contig_file, header = T, stringsAsFactors = F)
isos <- annots[,c('contig_id', 'c_gene')]
colnames(isos) <- c('ReadID_Heavy', 'Isotype')
umis <- annots[,c('contig_id', 'umis')]

### Heavy

heavy <- read.csv("Heavy/VH_RecombinationSummaries.RF.clones.txt", sep = '\t', header = T, stringsAsFactors = F)
heavy$ReadID <- heavy$UID
names(heavy)[names(heavy) == "Clone"] <- "CloneID"
heavy$ID <- gsub('_contig_[0-9]{1,2}', '', heavy$ReadID)
heavy$Heavy_Fxnl <- heavy$ReadID %in% heavy_fxnl$UID

heavy <- heavy %>%
  group_by(CloneID) %>%
  mutate(X.Members = n()) %>%
  ungroup()

# Rename all columns except 'ID' in the heavy dataframe
colnames(heavy)[colnames(heavy) != "ID"] <- paste0(colnames(heavy)[colnames(heavy) != "ID"], "_Heavy")
all_merged <- merge(heavy, chain_counts, by = 'ID', all.x = T)
all_merged <- merge(all_merged, isos, by = 'ReadID_Heavy', all.x = T)
all_merged <- merge(all_merged, umis, by.x = 'ReadID_Heavy', by.y = 'contig_id', all.x = T)
names(all_merged)[names(all_merged) == "umis"] <- "umis_Heavy"
names(all_merged)[names(all_merged) == "Heavy_Fxnl_Heavy"] <- "Heavy_Fxnl"
all_merged <- all_merged[order(all_merged$CloneID_Heavy),]

if (sum(is.na(all_merged$CDR3_Heavy))>=1){
  no_na <- subset(all_merged, !is.na(CDR3_Heavy))
  na <- subset(all_merged, is.na(CDR3_Heavy))
  no_na$CDR3_Heavy_AA <- Biostrings::translate(DNAStringSet(no_na$CDR3_Heavy))
  na$CDR3_Heavy_AA <- NA
  all_merged <- smartbind(no_na, na)
}else{
  all_merged$CDR3_Heavy_AA <- Biostrings::translate(DNAStringSet(all_merged$CDR3_Heavy))
}


all_merged$VGene_Heavy <- gsub('\\*[0-9]{1,}', '', all_merged$VGene_Heavy)
all_merged$VGene_Heavy <- gsub('IGHV', 'VH', all_merged$VGene_Heavy)
all_merged$DGene_Heavy <- gsub('\\*[0-9]{1,}', '', all_merged$DGene_Heavy)
all_merged$DGene_Heavy <- gsub('IGHD', 'DH', all_merged$DGene_Heavy)
all_merged$JGene_Heavy <- gsub('\\*[0-9]{1,}', '', all_merged$JGene_Heavy)
all_merged$JGene_Heavy <- gsub('IGHJ', 'JH', all_merged$JGene_Heavy)

all_merged <- select(all_merged, ID, CloneID_Heavy, X.Members_Heavy, ReadID_Heavy, VGene_Heavy, DGene_Heavy, JGene_Heavy, CDR3Length_Heavy, Isotype, MuFreq_Heavy, Heavy_Fxnl, N_Heavy, umis_Heavy, CDR3_Heavy, CDR3_Heavy_AA)

write.csv(all_merged, file = "10x_merged_clones.csv", row.names = F, quote = F)
