library(Seurat)
library(SingleR)
library(SeuratObject)
library(ggplot2)

setwd('/datacommons/dhvi/10X_Genomics_data/SHIV_BG505.T33N_Infection/Reverse_BEAM/CR7.2/Analysis/GEX/Seurat/1_preprocessing/')

# read in data and create seurat object
memoryPBMC_data <- Read10X('/datacommons/dhvi/10X_Genomics_data/SHIV_BG505.T33N_Infection/Reverse_BEAM/CR7.2/RM08N021/wk105/memory_PBMC/GEX/RM08N021_MemPBMC_Agg/outs/count/filtered_feature_bc_matrix/')
memoryPBMC <- CreateSeuratObject(counts = memoryPBMC_data, project = "memoryPBMC")

# basic qc
#memoryPBMC[["percent.mt"]] <- PercentageFeatureSet(memoryPBMC, pattern = "^MT-")
VlnPlot(memoryPBMC, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, alpha = .1)

memoryPBMC <- subset(memoryPBMC, subset = nFeature_RNA > 200 & nFeature_RNA < 6000)

# normalize, scale data
memoryPBMC <- NormalizeData(memoryPBMC)
memoryPBMC <- FindVariableFeatures(memoryPBMC, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(memoryPBMC)
memoryPBMC <- ScaleData(memoryPBMC, features = all.genes)

# read in fxnl 1:1 vdj data with BEAM data
vdj <- read.csv('/datacommons/dhvi/10X_Genomics_data/SHIV_BG505.T33N_Infection/Reverse_BEAM/CR7.2/RM08N021/wk105/memory_PBMC/BEAM/ShawHIVRAD_RM08N021wk105MemPBMC_merged_clonal_analysis.fxnl_1to1.BEAM_Scores.UMIs.Seqs.csv')
vdj <- subset(vdj, !(is.na(vdj$BG505_Mut11b_SOSIP_1_BEAM_Score)))
memoryPBMC@meta.data$VDJ_data <- ifelse(rownames(memoryPBMC@meta.data) %in% vdj$barcode, T, F)
table(memoryPBMC@meta.data$VDJ_data)

# run singleR to be able to subset only B cells
ref <- celldex::MonacoImmuneData()
bcell_subsets <- c('Exhausted B cells', 'Naive B cells', 'Non-switched memory B cells', 'Switched memory B cells', 'Plasmablasts')
singler.labels.main <- SingleR(as.SingleCellExperiment(memoryPBMC), ref = ref, labels = ref$label.main)
singler.labels.fine <- SingleR(as.SingleCellExperiment(memoryPBMC), ref = ref, labels = ref$label.fine)
memoryPBMC[["SingleR.labels.main"]] <- singler.labels.main$labels
memoryPBMC[["SingleR.labels.fine"]] <- singler.labels.fine$labels

# add vdj metadata
rownames(vdj) <- vdj$barcode
vdj <- dplyr::select(vdj, -barcode)
memoryPBMC <- AddMetaData(memoryPBMC, vdj)

# subset for only B cells (main and fine labels) with fxnl 1:1 vdj data
memoryPBMC_sub <- subset(memoryPBMC, SingleR.labels.main == 'B cells' & VDJ_data == T & SingleR.labels.fine %in% bcell_subsets)

# normalize, scale for integration with all samples
memoryPBMC_sub <- NormalizeData(memoryPBMC_sub)
memoryPBMC_sub <- FindVariableFeatures(memoryPBMC_sub, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(memoryPBMC_sub)
memoryPBMC_sub <- ScaleData(memoryPBMC_sub, features = all.genes)

save(memoryPBMC_sub, file = 'memoryPBMC_Bcells_wFxnl1to1VDJandBEAM_SeuratObject.RData')
nrow(memoryPBMC_sub@meta.data)

# memoryPBMC_sub <- subset(memoryPBMC, seurat_clusters != 3 & SingleR.labels.main == 'B cells')


