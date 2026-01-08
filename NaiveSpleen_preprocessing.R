library(Seurat)
library(SingleR)
library(SeuratObject)
library(ggplot2)

setwd('/datacommons/dhvi/10X_Genomics_data/SHIV_BG505.T33N_Infection/Reverse_BEAM/CR7.2/Analysis/GEX/Seurat/1_preprocessing/')

# read in data and create seurat object
naiveSpleen_data <- Read10X('/datacommons/dhvi/10X_Genomics_data/SHIV_BG505.T33N_Infection/Reverse_BEAM/CR7.2/RM08N021/wk105/naive_spleen/GEX/cellranger/Naive_Spleen_Agg/outs/count/filtered_feature_bc_matrix/')
naiveSpleen <- CreateSeuratObject(counts = naiveSpleen_data, project = "naiveSpleen")

# basic qc
#naiveSpleen[["percent.mt"]] <- PercentageFeatureSet(naiveSpleen, pattern = "^MT-")
VlnPlot(naiveSpleen, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, alpha = .1)

naiveSpleen <- subset(naiveSpleen, subset = nFeature_RNA > 200 & nFeature_RNA < 6000)

# normalize, scale data
naiveSpleen <- NormalizeData(naiveSpleen)
naiveSpleen <- FindVariableFeatures(naiveSpleen, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(naiveSpleen)
naiveSpleen <- ScaleData(naiveSpleen, features = all.genes)

# read in fxnl 1:1 vdj data with BEAM data
vdj <- read.csv('/datacommons/dhvi/10X_Genomics_data/SHIV_BG505.T33N_Infection/Reverse_BEAM/CR7.2/RM08N021/wk105/naive_spleen/BEAM/ShawHIVRAD_RM08N021wk105NaiveSpleen_merged_clonal_analysis.fxnl_1to1.BEAM_Scores.UMIs.Seqs.csv')
vdj <- subset(vdj, !(is.na(vdj$BG505_Mut11b_SOSIP_1_BEAM_Score)))
naiveSpleen@meta.data$VDJ_data <- ifelse(rownames(naiveSpleen@meta.data) %in% vdj$barcode, T, F)
table(naiveSpleen@meta.data$VDJ_data)

# run singleR to be able to subset only B cells
ref <- celldex::MonacoImmuneData()
bcell_subsets <- c('Exhausted B cells', 'Naive B cells', 'Non-switched memory B cells', 'Switched memory B cells', 'Plasmablasts')
singler.labels.main <- SingleR(as.SingleCellExperiment(naiveSpleen), ref = ref, labels = ref$label.main)
singler.labels.fine <- SingleR(as.SingleCellExperiment(naiveSpleen), ref = ref, labels = ref$label.fine)
naiveSpleen[["SingleR.labels.main"]] <- singler.labels.main$labels
naiveSpleen[["SingleR.labels.fine"]] <- singler.labels.fine$labels

# add vdj metadata
rownames(vdj) <- vdj$barcode
vdj <- dplyr::select(vdj, -barcode)
naiveSpleen <- AddMetaData(naiveSpleen, vdj)

# subset for only B cells (main and fine labels) with fxnl 1:1 vdj data
naiveSpleen_sub <- subset(naiveSpleen, SingleR.labels.main == 'B cells' & VDJ_data == T & SingleR.labels.fine %in% bcell_subsets)

# normalize, scale for integration with all samples
naiveSpleen_sub <- NormalizeData(naiveSpleen_sub)
naiveSpleen_sub <- FindVariableFeatures(naiveSpleen_sub, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(naiveSpleen_sub)
naiveSpleen_sub <- ScaleData(naiveSpleen_sub, features = all.genes)

save(naiveSpleen_sub, file = 'naiveSpleen_Bcells_wFxnl1to1VDJandBEAM_SeuratObject.RData')
nrow(naiveSpleen_sub@meta.data)

# naiveSpleen_sub <- subset(naiveSpleen, seurat_clusters != 3 & SingleR.labels.main == 'B cells')


