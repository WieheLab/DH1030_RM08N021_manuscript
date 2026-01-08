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
vdj_new <- vdj
vdj_new$barcode2 <- gsub('-[1-8]', '', vdj_new$barcode)
vdj_new$lane <- gsub('[AGCT]{16}-', '', vdj_new$barcode)

vdj_new$NewLane <- case_when(vdj_new$lane == '1' ~ '3',
                             vdj_new$lane == '2' ~ '2',
                             vdj_new$lane == '3' ~ '1',
                             vdj_new$lane == '4' ~ '4',
                             vdj_new$lane == '5' ~ '6',
                             vdj_new$lane == '6' ~ '5',
                             vdj_new$lane == '7' ~ '8',
                             vdj_new$lane == '8' ~ '7')

vdj_new$NewBarcode <- paste(vdj_new$barcode2, as.character(vdj_new$NewLane), sep = '-')
naiveLN@meta.data$VDJ_data <- ifelse(rownames(naiveLN@meta.data) %in% vdj_new$NewBarcode, T, F)
table(naiveLN@meta.data$VDJ_data)

# run singleR to be able to subset only B cells
ref <- celldex::MonacoImmuneData()
bcell_subsets <- c('Exhausted B cells', 'Naive B cells', 'Non-switched memory B cells', 'Switched memory B cells', 'Plasmablasts')
singler.labels.main <- SingleR(as.SingleCellExperiment(naiveLN), ref = ref, labels = ref$label.main)
singler.labels.fine <- SingleR(as.SingleCellExperiment(naiveLN), ref = ref, labels = ref$label.fine)
naiveLN[["SingleR.labels.main"]] <- singler.labels.main$labels
naiveLN[["SingleR.labels.fine"]] <- singler.labels.fine$labels

# add vdj metadata
rownames(vdj_new) <- vdj_new$NewBarcode
vdj_new <- dplyr::select(vdj_new, -NewBarcode, -NewLane, -barcode, -barcode2, -lane)
naiveLN <- AddMetaData(naiveLN, vdj_new)

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


