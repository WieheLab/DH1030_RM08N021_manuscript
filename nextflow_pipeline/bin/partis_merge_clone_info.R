.libPaths('/hpc/group/dhvi/WieheLab/mb488/R_Libs/R_v4.4.3')

library(stringr)
library(dplyr)
library(optparse)
library(gtools)
library(Biostrings)
library(seqinr)

args <- commandArgs(trailingOnly = TRUE)

contig_file <- args[1]

chain_counts <- read.csv("../intermediate_files/chain_counts.csv")
heavy_fxnl <- read.csv("../intermediate_files/heavies.fxnl.ids.txt", header = T, stringsAsFactors = F)
kappa_fxnl <- read.csv("../intermediate_files/kappas.fxnl.ids.txt", header = T, stringsAsFactors = F)
lambda_fxnl <- read.csv("../intermediate_files/lambdas.fxnl.ids.txt", header = T, stringsAsFactors = F)
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

heavy <- heavy[,c("ID", "ReadID", "CloneID", "X.Members", "Heavy_Fxnl", "VGene", "DGene", "JGene", "CDR3", "CDR3Length", "MuFreq")]

### Kappa
kappa <- read.csv("Kappa/VK_RecombinationSummaries.RF.clones.txt", sep = '\t', header = T, stringsAsFactors = F)
kappa$ReadID <- kappa$UID
names(kappa)[names(kappa) == "Clone"] <- "CloneID"
kappa$ID <- gsub('_contig_[0-9]{1,2}', '', kappa$ReadID)
kappa$Light_Fxnl <- kappa$ReadID %in% kappa_fxnl$UID

kappa <- kappa %>%
  group_by(CloneID) %>%
  mutate(X.Members = n()) %>%
  ungroup()

kappa$Chain <- 'Kappa'

kappa <- kappa[,c("ID", "ReadID", "CloneID", "X.Members", "Chain", "Light_Fxnl", "VGene", "JGene", "CDR3", "CDR3Length", "MuFreq")]

kappa <- data.frame(kappa)

### Lambda
lambda <- read.csv("Lambda/VL_RecombinationSummaries.RF.clones.txt", sep = '\t', header = T, stringsAsFactors = F)
lambda$ReadID <- lambda$UID
names(lambda)[names(lambda) == "Clone"] <- "CloneID"
lambda$ID <- gsub('_contig_[0-9]{1,2}', '', lambda$ReadID)
lambda$Light_Fxnl <- lambda$ReadID %in% lambda_fxnl$UID

lambda <- lambda %>%
  group_by(CloneID) %>%
  mutate(X.Members = n()) %>%
  ungroup()

lambda$Chain <- 'Lambda'

lambda <- lambda[,c("ID","ReadID", "CloneID", "X.Members", "Chain", "Light_Fxnl", "VGene", "JGene", "CDR3", "CDR3Length", "MuFreq")]

lambda <- data.frame(lambda)

### Lights 
lights <- smartbind(kappa, lambda)

all_merged <- merge(heavy, lights, by = 'ID', suffixes = c('_Heavy', '_Light'), all = T)
all_merged <- merge(all_merged, chain_counts, by = 'ID', all.x = T)
all_merged <- merge(all_merged, isos, by = 'ReadID_Heavy', all.x = T)
all_merged <- merge(all_merged, umis, by.x = 'ReadID_Heavy', by.y = 'contig_id', all.x = T)
all_merged <- merge(all_merged, umis, by.x = 'ReadID_Light', by.y = 'contig_id', all.x = T, suffixes = c('_Heavy', '_Light'))
all_merged <- all_merged[order(all_merged$CloneID_Heavy, all_merged$CloneID_Light),]


if (sum(is.na(all_merged$CDR3_Heavy))>=1){
  no_na <- subset(all_merged, !is.na(CDR3_Heavy))
  na <- subset(all_merged, is.na(CDR3_Heavy))
  no_na$CDR3_Heavy_AA <- Biostrings::translate(DNAStringSet(no_na$CDR3_Heavy))
  na$CDR3_Heavy_AA <- NA
  all_merged <- smartbind(no_na, na)
}else{
  all_merged$CDR3_Heavy_AA <- Biostrings::translate(DNAStringSet(all_merged$CDR3_Heavy))
}


if (sum(is.na(all_merged$CDR3_Light))>=1){
  no_na <- subset(all_merged, !is.na(CDR3_Light))
  na <- subset(all_merged, is.na(CDR3_Light))
  no_na$CDR3_Light_AA <- Biostrings::translate(DNAStringSet(no_na$CDR3_Light))
  na$CDR3_Light_AA <- NA
  all_merged <- smartbind(no_na, na)
}else{
all_merged$CDR3_Light_AA <- Biostrings::translate(DNAStringSet(all_merged$CDR3_Light))
}

all_merged$VGene_Heavy <- gsub('\\*[0-9]{1,}', '', all_merged$VGene_Heavy)
all_merged$VGene_Heavy <- gsub('IGHV', 'VH', all_merged$VGene_Heavy)
all_merged$DGene <- gsub('\\*[0-9]{1,}', '', all_merged$DGene)
all_merged$DGene <- gsub('IGHD', 'DH', all_merged$DGene)
all_merged$JGene_Heavy <- gsub('\\*[0-9]{1,}', '', all_merged$JGene_Heavy)
all_merged$JGene_Heavy <- gsub('IGHJ', 'JH', all_merged$JGene_Heavy)

all_merged$VGene_Light <- gsub('\\*[0-9]{1,}', '', all_merged$VGene_Light)
all_merged$VGene_Light <- gsub('IGLV', 'VL', all_merged$VGene_Light)
all_merged$VGene_Light <- gsub('IGKV', 'VK', all_merged$VGene_Light)
all_merged$JGene_Light <- gsub('\\*[0-9]{1,}', '', all_merged$JGene_Light)
all_merged$JGene_Light <- gsub('IGLJ', 'JL', all_merged$JGene_Light)
all_merged$JGene_Light <- gsub('IGKJ', 'JK', all_merged$JGene_Light)
all_merged$Fxnl_1to1 <- ifelse(all_merged$N_Heavy == 1 & all_merged$N_Light == 1 & all_merged$Heavy_Fxnl == 1 & all_merged$Light_Fxnl == 1, T, F)

# Add VDJ sequences 

heavy_seqs <- read.fasta('../unfiltered/Heavy/heavies.VDJ.fasta', as.string = T)
kappa_seqs <- read.fasta('../unfiltered/Kappa/kappas.VJ.fasta', as.string = T)
lambda_seqs <- read.fasta('../unfiltered/Lambda/lambdas.VJ.fasta', as.string = T)

heavies_df <- data.frame(ReadID_Heavy = names(heavy_seqs), Heavy_VDJ_Seq = toupper(as.character(heavy_seqs)))
kappas_df <- data.frame(ReadID_Light = names(kappa_seqs), Light_VJ_Seq = toupper(as.character(kappa_seqs)))
lambdas_df <- data.frame(ReadID_Light = names(lambda_seqs), Light_VJ_Seq = toupper(as.character(lambda_seqs)))

lights_df <- smartbind(kappas_df, lambdas_df)

all_merged <- merge(all_merged, heavies_df, by.x = 'ReadID_Heavy', all.x = T)
all_merged <- merge(all_merged, lights_df, by.x = 'ReadID_Light', all.x = T)

# Barcode hashing
#converted_barcodes <- read.csv('/hpc/group/dhvi/WieheLab/mb488/scripts/10x/converted_10x_barcodes.csv', header = TRUE)

#all_converted <- all_merged

#all_converted$Gem_Well <- gsub('[AGCT]{16}-', '', all_merged$ID)
#all_converted$Light_chain_convert <- case_when(
  #all_converted$Chain == 'Kappa' ~ 'K',
  #all_converted$Chain == 'Lambda' ~ 'L'
#)
#all_converted$Heavy_contig <- gsub('[ACGT]{16}-[1-8]_contig_', '', all_converted$ReadID_Heavy)
#all_converted$Light_contig <- gsub('[ACGT]{16}-[1-8]_contig_', '', all_converted$ReadID_Light)

#all_converted <- merge(all_converted, converted_barcodes, by.x = 'ID', by.y = 'X10X_BC')

#all_converted$ConvertedID <- paste(Exp, PTID, paste('G', all_converted$Gem_Well, sep = ''), all_converted$Converted_BC, sep = '_')
#all_converted$Converted_ReadID_Heavy <- paste(Exp, PTID, paste('G', all_converted$Gem_Well, sep = ''), all_converted$Converted_BC, paste('H', all_converted$Heavy_contig, sep = ''), sep = '_')
#all_converted$Converted_ReadID_Light <- paste(Exp, PTID, paste('G', all_converted$Gem_Well, sep = ''), all_converted$Converted_BC, paste(all_converted$Light_chain_convert, all_converted$Light_contig, sep = ''), sep = '_')


#all_converted <- select(all_converted, ID, Converted_ID, Fxnl_1to1, CloneID_Heavy, X.Members_Heavy, ReadID_Heavy, Converted_ReadID_Heavy, VGene_Heavy, DGene, JGene_Heavy, CDR3Length_Heavy, Isotype, MuFreq_Heavy,
                     #CloneID_Light, X.Members_Light, ReadID_Light, Converted_ReadID_Light, Chain, VGene_Light, JGene_Light, CDR3Length_Light, MuFreq_Light,
                     #Heavy_Fxnl, Light_Fxnl, N_Heavy, N_Kappa, N_Lambda, N_Light, umis_Heavy, umis_Light, CDR3_Heavy, CDR3_Heavy_AA, CDR3_Light, CDR3_Light_AA, Heavy_VDJ_Seq, Light_VJ_Seq)


#all_converted_subset <- subset(all_converted, N_Heavy == 1 & N_Light == 1 & Heavy_Fxnl == 1 & Light_Fxnl == 1)

all_merged <- select(all_merged, ID, Fxnl_1to1, CloneID_Heavy, X.Members_Heavy, ReadID_Heavy, VGene_Heavy, DGene, JGene_Heavy, CDR3Length_Heavy, Isotype, MuFreq_Heavy,
                     CloneID_Light, X.Members_Light, ReadID_Light, Chain, VGene_Light, JGene_Light, CDR3Length_Light, MuFreq_Light,
                     Heavy_Fxnl, Light_Fxnl, N_Heavy, N_Kappa, N_Lambda, N_Light, umis_Heavy, umis_Light, CDR3_Heavy, CDR3_Heavy_AA, CDR3_Light, CDR3_Light_AA, Heavy_VDJ_Seq, Light_VJ_Seq)
                     
all_merged_subset <- subset(all_merged, N_Heavy == 1 & N_Light == 1 & Heavy_Fxnl == 1 & Light_Fxnl == 1)

write.csv(all_merged, file = "10x_merged_clones.csv", row.names = F, quote = F)
write.csv(all_merged_subset, file = "10x_merged_clones.fxnl_1to1.csv", row.names = F, quote = F)