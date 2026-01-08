library(dplyr)
library(gtools)
library(seqinr)

args <- commandArgs(trailingOnly = TRUE)

VDJ_PATH <- args[1]
PARTIS_PATH <- args[2]

setwd(VDJ_PATH)

vdj <- read.csv(paste0(VDJ_PATH, "/clonal_analysis/10x_merged_clones.csv"))

partis <- read.csv(paste0(PARTIS_PATH, "/analysis/clonal_analysis/VH_RecombinationSummaries.txt"), sep='\t', header = T, stringsAsFactors = F)

partis_genes <- select(partis, UID, VGene_Heavy, DGene_Heavy, JGene_Heavy)
colnames(partis_genes) <- c('ReadID_Heavy', 'KimDB_Partis_VGene', 'KimDB_Partis_DGene', 'KimDB_Partis_JGene')

vdj_with_partis <- merge(vdj, partis_genes, by = 'ReadID_Heavy', all.x = T)
vdj_with_partis <- select(vdj_with_partis, ID:X.Members_Heavy, ReadID_Heavy, VGene_Heavy, KimDB_Partis_VGene, DGene, 
                          KimDB_Partis_DGene, JGene_Heavy, KimDB_Partis_JGene, CDR3Length_Heavy:CDR3_Light_AA)

heavies <- read.fasta('unfiltered/Heavy/heavies.VDJ.fasta', as.string = T, forceDNAtolower = F)
heavies_df <- data.frame(ReadID_Heavy = names(heavies), Heavy_VDJ_Seq = as.character(heavies))

kappas <- read.fasta('unfiltered/Kappa/kappas.VJ.fasta', as.string = T, forceDNAtolower = F)
kappas_df <- data.frame(ReadID_Light = names(kappas), Light_VJ_Seq = as.character(kappas))

lambdas <- read.fasta('unfiltered/Lambda/lambdas.VJ.fasta', as.string = T, forceDNAtolower = F)
lambdas_df <- data.frame(ReadID_Light = names(lambdas), Light_VJ_Seq = as.character(lambdas))

lights_df <- smartbind(kappas_df, lambdas_df)

vdj_with_partis <- merge(vdj_with_partis, heavies_df, by = 'ReadID_Heavy', all.x = T)
vdj_with_partis <- merge(vdj_with_partis, lights_df, by = 'ReadID_Light', all.x = T)

vdj_with_partis <- select(vdj_with_partis, ID:X.Members_Heavy, ReadID_Heavy, VGene_Heavy:X.Members_Light, ReadID_Light, Chain:Light_VJ_Seq)
vdj_with_partis <- vdj_with_partis[order(vdj_with_partis$CloneID_Heavy, vdj_with_partis$CloneID_Light),]

colnames(vdj_with_partis)[which(colnames(vdj_with_partis)=='DGene')] <- 'DGene_Heavy'

write.csv(vdj_with_partis, file = '../merged_clonal_analysis_with_KimDB_partis.csv', row.names = F, quote = F)
