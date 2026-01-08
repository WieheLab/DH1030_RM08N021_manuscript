library(Seurat)
library(SingleR)
library(SeuratObject)
library(ggplot2)

setwd('/datacommons/dhvi/10X_Genomics_data/SHIV_BG505.T33N_Infection/Reverse_BEAM/CR7.2/Analysis/GEX/Seurat/1_preprocessing/')

# read in data and create seurat object
gcLN_data <- Read10X('/datacommons/dhvi/10X_Genomics_data/SHIV_BG505.T33N_Infection/Reverse_BEAM/CR7.2/RM08N021/wk105/gc_LN/GEX/RM08N021_gcLN_Agg/outs/count/filtered_feature_bc_matrix/')
gcLN <- CreateSeuratObject(counts = gcLN_data, project = "gcLN")

# basic qc
#gcLN[["percent.mt"]] <- PercentageFeatureSet(gcLN, pattern = "^MT-")
VlnPlot(gcLN, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, alpha = .1)

gcLN <- subset(gcLN, subset = nFeature_RNA > 200 & nFeature_RNA < 6000)

# normalize, scale data
gcLN <- NormalizeData(gcLN)
gcLN <- FindVariableFeatures(gcLN, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(gcLN)
gcLN <- ScaleData(gcLN, features = all.genes)

# read in fxnl 1:1 vdj data with BEAM data
vdj <- read.csv('/datacommons/dhvi/10X_Genomics_data/SHIV_BG505.T33N_Infection/Reverse_BEAM/CR7.2/RM08N021/wk105/gc_LN/BEAM/ShawHIVRAD_RM08N021wk105GCLN_merged_clonal_analysis.fxnl_1to1.BEAM_Scores.UMIs.Seqs.csv')
vdj <- subset(vdj, !(is.na(vdj$BG505_Mut11b_SOSIP_1_BEAM_Score)))
gcLN@meta.data$VDJ_data <- ifelse(rownames(gcLN@meta.data) %in% vdj$barcode, T, F)
table(gcLN@meta.data$VDJ_data)

# run singleR to be able to subset only B cells
ref <- celldex::MonacoImmuneData()
bcell_subsets <- c('Exhausted B cells', 'Naive B cells', 'Non-switched memory B cells', 'Switched memory B cells', 'Plasmablasts')
singler.labels.main <- SingleR(as.SingleCellExperiment(gcLN), ref = ref, labels = ref$label.main)
singler.labels.fine <- SingleR(as.SingleCellExperiment(gcLN), ref = ref, labels = ref$label.fine)
gcLN[["SingleR.labels.main"]] <- singler.labels.main$labels
gcLN[["SingleR.labels.fine"]] <- singler.labels.fine$labels

# add vdj metadata
rownames(vdj) <- vdj$barcode
vdj <- dplyr::select(vdj, -barcode)
gcLN <- AddMetaData(gcLN, vdj)

# subset for only B cells (main and fine labels) with fxnl 1:1 vdj data
gcLN_sub <- subset(gcLN, SingleR.labels.main == 'B cells' & VDJ_data == T & SingleR.labels.fine %in% bcell_subsets)

# normalize, scale for integration with all samples
gcLN_sub <- NormalizeData(gcLN_sub)
gcLN_sub <- FindVariableFeatures(gcLN_sub, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(gcLN_sub)
gcLN_sub <- ScaleData(gcLN_sub, features = all.genes)

save(gcLN_sub, file = 'gcLN_Bcells_wFxnl1to1VDJandBEAM_SeuratObject.RData')
nrow(gcLN_sub@meta.data)

# gcLN_sub <- subset(gcLN, seurat_clusters != 3 & SingleR.labels.main == 'B cells')


