library(Seurat)
library(SingleR)
library(SeuratObject)
library(ggplot2)

setwd('/datacommons/dhvi/10X_Genomics_data/SHIV_BG505.T33N_Infection/Reverse_BEAM/CR7.2/Analysis/GEX/Seurat/1_preprocessing/')

# read in data and create seurat object
naiveLN_data <- Read10X('/datacommons/dhvi/10X_Genomics_data/SHIV_BG505.T33N_Infection/Reverse_BEAM/CR7.2/RM08N021/wk105/naive_LN/GEX/RM08N021_NaiveLN_Agg/outs/count/filtered_feature_bc_matrix/')
naiveLN <- CreateSeuratObject(counts = naiveLN_data, project = "naiveLN")

# basic qc
#naiveLN[["percent.mt"]] <- PercentageFeatureSet(naiveLN, pattern = "^MT-")
VlnPlot(naiveLN, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, alpha = .1)

naiveLN <- subset(naiveLN, subset = nFeature_RNA > 200 & nFeature_RNA < 6000)

# normalize, scale data
naiveLN <- NormalizeData(naiveLN)
naiveLN <- FindVariableFeatures(naiveLN, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(naiveLN)
naiveLN <- ScaleData(naiveLN, features = all.genes)

# read in fxnl 1:1 vdj data with BEAM data
vdj <- read.csv('/datacommons/dhvi/10X_Genomics_data/SHIV_BG505.T33N_Infection/Reverse_BEAM/CR7.2/RM08N021/wk105/naive_LN/BEAM/ShawHIVRAD_RM08N021wk105NaiveLN_merged_clonal_analysis.fxnl_1to1.BEAM_Scores.UMIs.Seqs.csv')
vdj <- subset(vdj, !(is.na(vdj$BG505_Mut11b_SOSIP_1_BEAM_Score)))
naiveLN@meta.data$VDJ_data <- ifelse(rownames(naiveLN@meta.data) %in% vdj$barcode, T, F)
table(naiveLN@meta.data$VDJ_data)

# run singleR to be able to subset only B cells
ref <- celldex::MonacoImmuneData()
bcell_subsets <- c('Exhausted B cells', 'Naive B cells', 'Non-switched naive B cells', 'Switched naive B cells', 'Plasmablasts')
singler.labels.main <- SingleR(as.SingleCellExperiment(naiveLN), ref = ref, labels = ref$label.main)
singler.labels.fine <- SingleR(as.SingleCellExperiment(naiveLN), ref = ref, labels = ref$label.fine)
naiveLN[["SingleR.labels.main"]] <- singler.labels.main$labels
naiveLN[["SingleR.labels.fine"]] <- singler.labels.fine$labels

# add vdj metadata
rownames(vdj) <- vdj$barcode
vdj <- dplyr::select(vdj, -barcode)
naiveLN <- AddMetaData(naiveLN, vdj)

# subset for only B cells (main and fine labels) with fxnl 1:1 vdj data
naiveLN_sub <- subset(naiveLN, SingleR.labels.main == 'B cells' & VDJ_data == T & SingleR.labels.fine %in% bcell_subsets)

# normalize, scale for integration with all samples
naiveLN_sub <- NormalizeData(naiveLN_sub)
naiveLN_sub <- FindVariableFeatures(naiveLN_sub, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(naiveLN_sub)
naiveLN_sub <- ScaleData(naiveLN_sub, features = all.genes)

save(naiveLN_sub, file = 'naiveLN_Bcells_wFxnl1to1VDJandBEAM_SeuratObject.RData')
nrow(naiveLN_sub@meta.data)

# naiveLN_sub <- subset(naiveLN, seurat_clusters != 3 & SingleR.labels.main == 'B cells')


