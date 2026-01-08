library(Seurat)
library(SingleR)
library(SeuratObject)
library(ggplot2)

setwd('/datacommons/dhvi/10X_Genomics_data/SHIV_BG505.T33N_Infection/Reverse_BEAM/CR7.2/Analysis/GEX/Seurat/1_preprocessing/')

# read in data and create seurat object
memoryLN_data <- Read10X('/datacommons/dhvi/10X_Genomics_data/SHIV_BG505.T33N_Infection/Reverse_BEAM/CR7.2/RM08N021/wk105/memory_LN/GEX/cellranger/RM08N021_MemLN_Agg/outs/count/filtered_feature_bc_matrix/')
memoryLN <- CreateSeuratObject(counts = memoryLN_data, project = "memoryLN")

# basic qc
#memoryLN[["percent.mt"]] <- PercentageFeatureSet(memoryLN, pattern = "^MT-")
VlnPlot(memoryLN, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, alpha = .1)

memoryLN <- subset(memoryLN, subset = nFeature_RNA > 200 & nFeature_RNA < 6000)

# normalize, scale data
memoryLN <- NormalizeData(memoryLN)
memoryLN <- FindVariableFeatures(memoryLN, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(memoryLN)
memoryLN <- ScaleData(memoryLN, features = all.genes)

# read in fxnl 1:1 vdj data with BEAM data
vdj <- read.csv('/datacommons/dhvi/10X_Genomics_data/SHIV_BG505.T33N_Infection/Reverse_BEAM/CR7.2/RM08N021/wk105/memory_LN/BEAM/ShawHIVRAD_RM08N021wk105MemLN_merged_clonal_analysis.fxnl_1to1.BEAM_Scores.UMIs.Seqs.csv')
vdj <- subset(vdj, !(is.na(vdj$BG505_Mut11b_SOSIP_1_BEAM_Score)))
memoryLN@meta.data$VDJ_data <- ifelse(rownames(memoryLN@meta.data) %in% vdj$barcode, T, F)
table(memoryLN@meta.data$VDJ_data)

# run singleR to be able to subset only B cells
ref <- celldex::MonacoImmuneData()
bcell_subsets <- c('Exhausted B cells', 'Naive B cells', 'Non-switched memory B cells', 'Switched memory B cells', 'Plasmablasts')
singler.labels.main <- SingleR(as.SingleCellExperiment(memoryLN), ref = ref, labels = ref$label.main)
singler.labels.fine <- SingleR(as.SingleCellExperiment(memoryLN), ref = ref, labels = ref$label.fine)
memoryLN[["SingleR.labels.main"]] <- singler.labels.main$labels
memoryLN[["SingleR.labels.fine"]] <- singler.labels.fine$labels

# add vdj metadata
rownames(vdj) <- vdj$barcode
vdj <- dplyr::select(vdj, -barcode)
memoryLN <- AddMetaData(memoryLN, vdj)

# subset for only B cells (main and fine labels) with fxnl 1:1 vdj data
memoryLN_sub <- subset(memoryLN, SingleR.labels.main == 'B cells' & VDJ_data == T & SingleR.labels.fine %in% bcell_subsets)

# normalize, scale for integration with all samples
memoryLN_sub <- NormalizeData(memoryLN_sub)
memoryLN_sub <- FindVariableFeatures(memoryLN_sub, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(memoryLN_sub)
memoryLN_sub <- ScaleData(memoryLN_sub, features = all.genes)

save(memoryLN_sub, file = 'memoryLN_Bcells_wFxnl1to1VDJandBEAM_SeuratObject.RData')
nrow(memoryLN_sub@meta.data)

# memoryLN_sub <- subset(memoryLN, seurat_clusters != 3 & SingleR.labels.main == 'B cells')


