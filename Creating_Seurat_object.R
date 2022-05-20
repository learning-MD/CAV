library(Seurat)
library(tidyverse)
library(patchwork)
library(DropletUtils)
library(glmGamPoi)
library(sctransform)
library(harmony)
library(SeuratDisk)
memory.limit(10000000000)

## Sample 1
ctrl_1 <- Read10X_h5("6460_KA_1/KA_1_output.h5")
ctrl_1 <- CreateSeuratObject(counts = ctrl_1, min.cells = 3, min.features = 100)
ctrl_1[["percent.mt"]] <- PercentageFeatureSet(ctrl_1, pattern = "^MT-")

# Subset features for [3rd quartile + 1.5(IQR)] the nCounts and nFeatures
quantile(ctrl_1$nFeature_RNA)
IQR(ctrl_1$nFeature_RNA)
1.5 * 869.5
1304.25 + 1342.75

quantile(ctrl_1$nCount_RNA)
IQR(ctrl_1$nCount_RNA)
1.5 * 1516
2274 + 2095

# Pre-processing of the Seurat object
ctrl_1 <- subset(ctrl_1, subset = nFeature_RNA > 100 & nFeature_RNA < 2647 & nCount_RNA < 4369 & percent.mt < 5)

# Write outputs for scrublet analysis
write10xCounts(ctrl_1@assays$RNA@counts, path = "ctrl_1.output")

# Integrate scrublet results back into SeuratObject
ctrl_1_scrublet <- read_tsv("scrublet_results/ctrl_1_scrublet_calls.tsv", col_names = FALSE)
ctrl_1_scrublet <- dplyr::rename(ctrl_1_scrublet, cell_ID = X1, doublet_score = X2, doublet = X3)
merged_metadata <- ctrl_1@meta.data
merged_metadata <- rownames_to_column(merged_metadata, var = "cell_ID")
merged_metadata <- left_join(merged_metadata, ctrl_1_scrublet, by = "cell_ID")
merged_metadata <- column_to_rownames(merged_metadata, var = "cell_ID")
ctrl_1@meta.data <- merged_metadata

# Souporcell
View(ctrl_1@meta.data)
ctrl_1_souporcell <- read_tsv("KA_1_souporcell_clusters.tsv")
ctrl_1@meta.data <- rownames_to_column(ctrl_1@meta.data, var = "barcode")
ctrl_1_merged_metadata <- ctrl_1@meta.data
ctrl_1_merged_metadata <- left_join(ctrl_1_merged_metadata, ctrl_1_souporcell, by = "barcode")
ctrl_1_merged_metadata <- column_to_rownames(ctrl_1_merged_metadata, var = "barcode")
ctrl_1@meta.data <- ctrl_1_merged_metadata


ctrl_1 <- subset(ctrl_1, subset = doublet == "FALSE")
ctrl_1@meta.data$cells <- "RV_control"
ctrl_1@meta.data$sample_ID <- "6460_KA_1"

## Sample 2

ctrl_2 <- Read10X_h5("6460_KA_2/KA_2_output.h5")
ctrl_2 <- CreateSeuratObject(counts = ctrl_2, min.cells = 3, min.features = 100)
ctrl_2[["percent.mt"]] <- PercentageFeatureSet(ctrl_2, pattern = "^MT-")

# Subset features for [3rd quartile + 1.5(IQR)] the nCounts and nFeatures
quantile(ctrl_2$nFeature_RNA)
IQR(ctrl_2$nFeature_RNA)
1.5 * 681
1021.5 + 1098

quantile(ctrl_2$nCount_RNA)
IQR(ctrl_2$nCount_RNA)
1.5 * 1049.75
1574.625 + 1551.75

# Pre-processing of the Seurat object
ctrl_2 <- subset(ctrl_2, subset = nFeature_RNA > 100 & nFeature_RNA < 2120 & nCount_RNA < 3127 & percent.mt < 5)

# Write outputs for scrublet analysis
write10xCounts(ctrl_2@assays$RNA@counts, path = "ctrl_2.output")

# Integrate scrublet results back into SeuratObject
ctrl_2_scrublet <- read_tsv("scrublet_results/ctrl_2_scrublet_calls.tsv", col_names = FALSE)
ctrl_2_scrublet <- dplyr::rename(ctrl_2_scrublet, cell_ID = X1, doublet_score = X2, doublet = X3)
ctrl_2_merged_metadata <- ctrl_2@meta.data
ctrl_2_merged_metadata <- rownames_to_column(ctrl_2_merged_metadata, var = "cell_ID")
ctrl_2_merged_metadata <- left_join(ctrl_2_merged_metadata, ctrl_2_scrublet, by = "cell_ID")
ctrl_2_merged_metadata <- column_to_rownames(ctrl_2_merged_metadata, var = "cell_ID")
ctrl_2@meta.data <- ctrl_2_merged_metadata

# Souporcell
ctrl_2_souporcell <- read_tsv("KA_2_souporcell_clusters.tsv")
ctrl_2@meta.data <- rownames_to_column(ctrl_2@meta.data, var = "barcode")
ctrl_2_merged_metadata <- ctrl_2@meta.data
ctrl_2_merged_metadata <- left_join(ctrl_2_merged_metadata, ctrl_2_souporcell, by = "barcode")
ctrl_2_merged_metadata <- column_to_rownames(ctrl_2_merged_metadata, var = "barcode")
ctrl_2@meta.data <- ctrl_2_merged_metadata


ctrl_2 <- subset(ctrl_2, subset = doublet == "FALSE")
ctrl_2@meta.data$cells <- "RV_control"
ctrl_2@meta.data$sample_ID <- "6460_KA_2"

## Sample 3

ctrl_3 <- Read10X_h5("6460_KA_3/KA_3_output.h5")
ctrl_3 <- CreateSeuratObject(counts = ctrl_3, min.cells = 3, min.features = 100)
ctrl_3[["percent.mt"]] <- PercentageFeatureSet(ctrl_3, pattern = "^MT-")

# Subset features for [3rd quartile + 1.5(IQR)] the nCounts and nFeatures
quantile(ctrl_3$nFeature_RNA)
IQR(ctrl_3$nFeature_RNA)
1.5 * 342
513 + 610

quantile(ctrl_3$nCount_RNA)
IQR(ctrl_3$nCount_RNA)
1.5 * 513.25
769.875 + 835.25

# Pre-processing of the Seurat object
ctrl_3 <- subset(ctrl_3, subset = nFeature_RNA > 100 & nFeature_RNA < 1123 & nCount_RNA < 1606 & percent.mt < 5)

# Write outputs for scrublet analysis
write10xCounts(ctrl_3@assays$RNA@counts, path = "ctrl_3.output")

# Integrate scrublet results back into SeuratObject
ctrl_3_scrublet <- read_tsv("scrublet_results/ctrl_3_scrublet_calls.tsv", col_names = FALSE)
ctrl_3_scrublet <- dplyr::rename(ctrl_3_scrublet, cell_ID = X1, doublet_score = X2, doublet = X3)
ctrl_3_merged_metadata <- ctrl_3@meta.data
ctrl_3_merged_metadata <- rownames_to_column(ctrl_3_merged_metadata, var = "cell_ID")
ctrl_3_merged_metadata <- left_join(ctrl_3_merged_metadata, ctrl_3_scrublet, by = "cell_ID")
ctrl_3_merged_metadata <- column_to_rownames(ctrl_3_merged_metadata, var = "cell_ID")
ctrl_3@meta.data <- ctrl_3_merged_metadata

# Souporcell
ctrl_3_souporcell <- read_tsv("KA_3_souporcell_clusters.tsv")
ctrl_3@meta.data <- rownames_to_column(ctrl_3@meta.data, var = "barcode")
ctrl_3_merged_metadata <- ctrl_3@meta.data
ctrl_3_merged_metadata <- left_join(ctrl_3_merged_metadata, ctrl_3_souporcell, by = "barcode")
ctrl_3_merged_metadata <- column_to_rownames(ctrl_3_merged_metadata, var = "barcode")
ctrl_3@meta.data <- ctrl_3_merged_metadata

ctrl_3 <- subset(ctrl_3, subset = doublet == "FALSE")
ctrl_3@meta.data$cells <- "RV_CAV"
ctrl_3@meta.data$sample_ID <- "6460_KA_3"

## Sample 4

ctrl_4 <- Read10X_h5("6460_KA_4/KA_4_output.h5")
ctrl_4 <- CreateSeuratObject(counts = ctrl_4, min.cells = 3, min.features = 100)
ctrl_4[["percent.mt"]] <- PercentageFeatureSet(ctrl_4, pattern = "^MT-")

# Subset features for [3rd quartile + 1.5(IQR)] the nCounts and nFeatures
quantile(ctrl_4$nFeature_RNA)
IQR(ctrl_4$nFeature_RNA)
1.5 * 1163
1744.5 + 1674

quantile(ctrl_4$nCount_RNA)
IQR(ctrl_4$nCount_RNA)
1.5 * 2504
3756 + 3163

# Pre-processing of the Seurat object
ctrl_4 <- subset(ctrl_4, subset = nFeature_RNA > 100 & nFeature_RNA < 3419 & nCount_RNA < 6919 & percent.mt < 5)

# Write outputs for scrublet analysis
write10xCounts(ctrl_4@assays$RNA@counts, path = "ctrl_4.output")

# Integrate scrublet results back into SeuratObject
ctrl_4_scrublet <- read_tsv("scrublet_results/ctrl_4_scrublet_calls.tsv", col_names = FALSE)
ctrl_4_scrublet <- dplyr::rename(ctrl_4_scrublet, cell_ID = X1, doublet_score = X2, doublet = X3)
ctrl_4_merged_metadata <- ctrl_4@meta.data
ctrl_4_merged_metadata <- rownames_to_column(ctrl_4_merged_metadata, var = "cell_ID")
ctrl_4_merged_metadata <- left_join(ctrl_4_merged_metadata, ctrl_4_scrublet, by = "cell_ID")
ctrl_4_merged_metadata <- column_to_rownames(ctrl_4_merged_metadata, var = "cell_ID")
ctrl_4@meta.data <- ctrl_4_merged_metadata

# Souporcell
ctrl_4_souporcell <- read_tsv("KA_4_souporcell_clusters.tsv")
ctrl_4@meta.data <- rownames_to_column(ctrl_4@meta.data, var = "barcode")
ctrl_4_merged_metadata <- ctrl_4@meta.data
ctrl_4_merged_metadata <- left_join(ctrl_4_merged_metadata, ctrl_4_souporcell, by = "barcode")
ctrl_4_merged_metadata <- column_to_rownames(ctrl_4_merged_metadata, var = "barcode")
ctrl_4@meta.data <- ctrl_4_merged_metadata

ctrl_4 <- subset(ctrl_4, subset = doublet == "FALSE")
ctrl_4@meta.data$cells <- "LV_CAV"
ctrl_4@meta.data$sample_ID <- "6460_KA_4"

## Sample 5

ctrl_5 <- Read10X_h5("6460_KA_5/KA_5_output.h5")
ctrl_5 <- CreateSeuratObject(counts = ctrl_5, min.cells = 3, min.features = 100)
ctrl_5[["percent.mt"]] <- PercentageFeatureSet(ctrl_5, pattern = "^MT-")

# Subset features for [3rd quartile + 1.5(IQR)] the nCounts and nFeatures
quantile(ctrl_5$nFeature_RNA)
IQR(ctrl_5$nFeature_RNA)
1.5 * 634
951 + 951

quantile(ctrl_5$nCount_RNA)
IQR(ctrl_5$nCount_RNA)
1.5 * 1029
1543.5 + 1402

# Pre-processing of the Seurat object
ctrl_5 <- subset(ctrl_5, subset = nFeature_RNA > 100 & nFeature_RNA < 1902 & nCount_RNA < 2946 & percent.mt < 5)

# Write outputs for scrublet analysis
write10xCounts(ctrl_5@assays$RNA@counts, path = "ctrl_5.output")

# Integrate scrublet results back into SeuratObject
ctrl_5_scrublet <- read_tsv("scrublet_results/ctrl_5_scrublet_calls.tsv", col_names = FALSE)
ctrl_5_scrublet <- dplyr::rename(ctrl_5_scrublet, cell_ID = X1, doublet_score = X2, doublet = X3)
ctrl_5_merged_metadata <- ctrl_5@meta.data
ctrl_5_merged_metadata <- rownames_to_column(ctrl_5_merged_metadata, var = "cell_ID")
ctrl_5_merged_metadata <- left_join(ctrl_5_merged_metadata, ctrl_5_scrublet, by = "cell_ID")
ctrl_5_merged_metadata <- column_to_rownames(ctrl_5_merged_metadata, var = "cell_ID")
ctrl_5@meta.data <- ctrl_5_merged_metadata

# Souporcell
ctrl_5_souporcell <- read_tsv("KA_5_souporcell_clusters.tsv")
ctrl_5@meta.data <- rownames_to_column(ctrl_5@meta.data, var = "barcode")
ctrl_5_merged_metadata <- ctrl_5@meta.data
ctrl_5_merged_metadata <- left_join(ctrl_5_merged_metadata, ctrl_5_souporcell, by = "barcode")
ctrl_5_merged_metadata <- column_to_rownames(ctrl_5_merged_metadata, var = "barcode")
ctrl_5@meta.data <- ctrl_5_merged_metadata

ctrl_5 <- subset(ctrl_5, subset = doublet == "FALSE")
ctrl_5@meta.data$cells <- "RV_control"
ctrl_5@meta.data$sample_ID <- "6460_KA_5"

## Sample 6

ctrl_6 <- Read10X_h5("6460_KA_6/KA_6_output.h5")
ctrl_6 <- CreateSeuratObject(counts = ctrl_6, min.cells = 3, min.features = 100)
ctrl_6[["percent.mt"]] <- PercentageFeatureSet(ctrl_6, pattern = "^MT-")

# Subset features for [3rd quartile + 1.5(IQR)] the nCounts and nFeatures
quantile(ctrl_6$nFeature_RNA)
IQR(ctrl_6$nFeature_RNA)
1.5 * 618
927 + 863

quantile(ctrl_6$nCount_RNA)
IQR(ctrl_6$nCount_RNA)
1.5 * 935
1402.5 + 1216

# Pre-processing of the Seurat object
ctrl_6 <- subset(ctrl_6, subset = nFeature_RNA > 100 & nFeature_RNA < 1790 & nCount_RNA < 2619 & percent.mt < 5)

# Write outputs for scrublet analysis
write10xCounts(ctrl_6@assays$RNA@counts, path = "ctrl_6.output")

# Integrate scrublet results back into SeuratObject
ctrl_6_scrublet <- read_tsv("scrublet_results/ctrl_6_scrublet_calls.tsv", col_names = FALSE)
ctrl_6_scrublet <- dplyr::rename(ctrl_6_scrublet, cell_ID = X1, doublet_score = X2, doublet = X3)
ctrl_6_merged_metadata <- ctrl_6@meta.data
ctrl_6_merged_metadata <- rownames_to_column(ctrl_6_merged_metadata, var = "cell_ID")
ctrl_6_merged_metadata <- left_join(ctrl_6_merged_metadata, ctrl_6_scrublet, by = "cell_ID")
ctrl_6_merged_metadata <- column_to_rownames(ctrl_6_merged_metadata, var = "cell_ID")
ctrl_6@meta.data <- ctrl_6_merged_metadata

# Souporcell
ctrl_6_souporcell <- read_tsv("KA_6_souporcell_clusters.tsv")
ctrl_6@meta.data <- rownames_to_column(ctrl_6@meta.data, var = "barcode")
ctrl_6_merged_metadata <- ctrl_6@meta.data
ctrl_6_merged_metadata <- left_join(ctrl_6_merged_metadata, ctrl_6_souporcell, by = "barcode")
ctrl_6_merged_metadata <- column_to_rownames(ctrl_6_merged_metadata, var = "barcode")
ctrl_6@meta.data <- ctrl_6_merged_metadata

ctrl_6 <- subset(ctrl_6, subset = doublet == "FALSE")
ctrl_6@meta.data$cells <- "RV_CAV"
ctrl_6@meta.data$sample_ID <- "6460_KA_6"

## Sample 7

ctrl_7 <- Read10X_h5("6460_KA_7/KA_7_output.h5")
ctrl_7 <- CreateSeuratObject(counts = ctrl_7, min.cells = 3, min.features = 100)
ctrl_7[["percent.mt"]] <- PercentageFeatureSet(ctrl_7, pattern = "^MT-")

# Subset features for [3rd quartile + 1.5(IQR)] the nCounts and nFeatures
quantile(ctrl_7$nFeature_RNA)
IQR(ctrl_7$nFeature_RNA)
1.5 * 587
880.5 + 942

quantile(ctrl_7$nCount_RNA)
IQR(ctrl_7$nCount_RNA)
1.5 * 945
1417.5 + 1376

# Pre-processing of the Seurat object
ctrl_7 <- subset(ctrl_7, subset = nFeature_RNA > 100 & nFeature_RNA < 1823 & nCount_RNA < 2794 & percent.mt < 5)

# Write outputs for scrublet analysis
write10xCounts(ctrl_7@assays$RNA@counts, path = "ctrl_7.output")

# Integrate scrublet results back into SeuratObject
ctrl_7_scrublet <- read_tsv("scrublet_results/ctrl_7_scrublet_calls.tsv", col_names = FALSE)
ctrl_7_scrublet <- dplyr::rename(ctrl_7_scrublet, cell_ID = X1, doublet_score = X2, doublet = X3)
ctrl_7_merged_metadata <- ctrl_7@meta.data
ctrl_7_merged_metadata <- rownames_to_column(ctrl_7_merged_metadata, var = "cell_ID")
ctrl_7_merged_metadata <- left_join(ctrl_7_merged_metadata, ctrl_7_scrublet, by = "cell_ID")
ctrl_7_merged_metadata <- column_to_rownames(ctrl_7_merged_metadata, var = "cell_ID")
ctrl_7@meta.data <- ctrl_7_merged_metadata

# Souporcell
ctrl_7_souporcell <- read_tsv("KA_7_souporcell_clusters.tsv")
ctrl_7@meta.data <- rownames_to_column(ctrl_7@meta.data, var = "barcode")
ctrl_7_merged_metadata <- ctrl_7@meta.data
ctrl_7_merged_metadata <- left_join(ctrl_7_merged_metadata, ctrl_7_souporcell, by = "barcode")
ctrl_7_merged_metadata <- column_to_rownames(ctrl_7_merged_metadata, var = "barcode")
ctrl_7@meta.data <- ctrl_7_merged_metadata

ctrl_7 <- subset(ctrl_7, subset = doublet == "FALSE")
ctrl_7@meta.data$cells <- "RV_CAV"
ctrl_7@meta.data$sample_ID <- "6460_KA_7"

## Sample 8

ctrl_8 <- Read10X_h5("6460_KA_8/KA_8_output.h5")
ctrl_8 <- CreateSeuratObject(counts = ctrl_8, min.cells = 3, min.features = 100)
ctrl_8[["percent.mt"]] <- PercentageFeatureSet(ctrl_8, pattern = "^MT-")

# Subset features for [3rd quartile + 1.5(IQR)] the nCounts and nFeatures
quantile(ctrl_8$nFeature_RNA)
IQR(ctrl_8$nFeature_RNA)
1.5 * 1553
2329.5 + 2575

quantile(ctrl_8$nCount_RNA)
IQR(ctrl_8$nCount_RNA)
1.5 * 4476.5
6714.75 + 6062.25

# Pre-processing of the Seurat object
ctrl_8 <- subset(ctrl_8, subset = nFeature_RNA > 100 & nFeature_RNA < 4905 & nCount_RNA < 12777 & percent.mt < 5)

# Write outputs for scrublet analysis
write10xCounts(ctrl_8@assays$RNA@counts, path = "ctrl_8.output")

# Integrate scrublet results back into SeuratObject
ctrl_8_scrublet <- read_tsv("scrublet_results/ctrl_8_scrublet_calls.tsv", col_names = FALSE)
ctrl_8_scrublet <- dplyr::rename(ctrl_8_scrublet, cell_ID = X1, doublet_score = X2, doublet = X3)
ctrl_8_merged_metadata <- ctrl_8@meta.data
ctrl_8_merged_metadata <- rownames_to_column(ctrl_8_merged_metadata, var = "cell_ID")
ctrl_8_merged_metadata <- left_join(ctrl_8_merged_metadata, ctrl_8_scrublet, by = "cell_ID")
ctrl_8_merged_metadata <- column_to_rownames(ctrl_8_merged_metadata, var = "cell_ID")
ctrl_8@meta.data <- ctrl_8_merged_metadata

# Souporcell
ctrl_8_souporcell <- read_tsv("KA_8_souporcell_clusters.tsv")
ctrl_8@meta.data <- rownames_to_column(ctrl_8@meta.data, var = "barcode")
ctrl_8_merged_metadata <- ctrl_8@meta.data
ctrl_8_merged_metadata <- left_join(ctrl_8_merged_metadata, ctrl_8_souporcell, by = "barcode")
ctrl_8_merged_metadata <- column_to_rownames(ctrl_8_merged_metadata, var = "barcode")
ctrl_8@meta.data <- ctrl_8_merged_metadata

ctrl_8 <- subset(ctrl_8, subset = doublet == "FALSE")
ctrl_8@meta.data$cells <- "LV_CAV"
ctrl_8@meta.data$sample_ID <- "6460_KA_8"

## Sample 9

ctrl_9 <- Read10X_h5("6460_KA_9/KA_9_output.h5")
ctrl_9 <- CreateSeuratObject(counts = ctrl_9, min.cells = 3, min.features = 100)
ctrl_9[["percent.mt"]] <- PercentageFeatureSet(ctrl_9, pattern = "^MT-")

# Subset features for [3rd quartile + 1.5(IQR)] the nCounts and nFeatures
quantile(ctrl_9$nFeature_RNA)
IQR(ctrl_9$nFeature_RNA)
1.5 * 589
883.5 + 996

quantile(ctrl_9$nCount_RNA)
IQR(ctrl_9$nCount_RNA)
1.5 * 932
1398 + 1434

# Pre-processing of the Seurat object
ctrl_9 <- subset(ctrl_9, subset = nFeature_RNA > 100 & nFeature_RNA < 1880 & nCount_RNA < 2832 & percent.mt < 5)

# Write outputs for scrublet analysis
write10xCounts(ctrl_9@assays$RNA@counts, path = "ctrl_9.output")

# Integrate scrublet results back into SeuratObject
ctrl_9_scrublet <- read_tsv("scrublet_results/ctrl_9_scrublet_calls.tsv", col_names = FALSE)
ctrl_9_scrublet <- dplyr::rename(ctrl_9_scrublet, cell_ID = X1, doublet_score = X2, doublet = X3)
ctrl_9_merged_metadata <- ctrl_9@meta.data
ctrl_9_merged_metadata <- rownames_to_column(ctrl_9_merged_metadata, var = "cell_ID")
ctrl_9_merged_metadata <- left_join(ctrl_9_merged_metadata, ctrl_9_scrublet, by = "cell_ID")
ctrl_9_merged_metadata <- column_to_rownames(ctrl_9_merged_metadata, var = "cell_ID")
ctrl_9@meta.data <- ctrl_9_merged_metadata

# Souporcell
ctrl_9_souporcell <- read_tsv("KA_9_souporcell_clusters.tsv")
ctrl_9_souporcell <- dplyr::rename(ctrl_9_souporcell, barcode = cell_ID)
ctrl_9@meta.data <- rownames_to_column(ctrl_9@meta.data, var = "barcode")
ctrl_9_merged_metadata <- ctrl_9@meta.data
ctrl_9_merged_metadata <- left_join(ctrl_9_merged_metadata, ctrl_9_souporcell, by = "barcode")
ctrl_9_merged_metadata <- column_to_rownames(ctrl_9_merged_metadata, var = "barcode")
ctrl_9@meta.data <- ctrl_9_merged_metadata

ctrl_9 <- subset(ctrl_9, subset = doublet == "FALSE")
ctrl_9@meta.data$cells <- "RV_CAV"
ctrl_9@meta.data$sample_ID <- "6460_KA_9"

## Sample 10

ctrl_10 <- Read10X_h5("6460_KA_10/KA_10_output.h5")
ctrl_10 <- CreateSeuratObject(counts = ctrl_10, min.cells = 3, min.features = 100)
ctrl_10[["percent.mt"]] <- PercentageFeatureSet(ctrl_10, pattern = "^MT-")

# Subset features for [3rd quartile + 1.5(IQR)] the nCounts and nFeatures
quantile(ctrl_10$nFeature_RNA)
IQR(ctrl_10$nFeature_RNA)
1.5 * 1157.25
1735.875 + 1592.25

quantile(ctrl_10$nCount_RNA)
IQR(ctrl_10$nCount_RNA)
1.5 * 2297
3445.5 + 2853

# Pre-processing of the Seurat object
ctrl_10 <- subset(ctrl_10, subset = nFeature_RNA > 100 & nFeature_RNA < 3329 & nCount_RNA < 6299 & percent.mt < 5)

# Write outputs for scrublet analysis
write10xCounts(ctrl_10@assays$RNA@counts, path = "ctrl_10.output")

# Integrate scrublet results back into SeuratObject
ctrl_10_scrublet <- read_tsv("scrublet_results/ctrl_10_scrublet_calls.tsv", col_names = FALSE)
ctrl_10_scrublet <- dplyr::rename(ctrl_10_scrublet, cell_ID = X1, doublet_score = X2, doublet = X3)
ctrl_10_merged_metadata <- ctrl_10@meta.data
ctrl_10_merged_metadata <- rownames_to_column(ctrl_10_merged_metadata, var = "cell_ID")
ctrl_10_merged_metadata <- left_join(ctrl_10_merged_metadata, ctrl_10_scrublet, by = "cell_ID")
ctrl_10_merged_metadata <- column_to_rownames(ctrl_10_merged_metadata, var = "cell_ID")
ctrl_10@meta.data <- ctrl_10_merged_metadata

# Souporcell
ctrl_10_souporcell <- read_tsv("KA_10_souporcell_clusters.tsv")
ctrl_10@meta.data <- rownames_to_column(ctrl_10@meta.data, var = "barcode")
ctrl_10_merged_metadata <- ctrl_10@meta.data
ctrl_10_merged_metadata <- left_join(ctrl_10_merged_metadata, ctrl_10_souporcell, by = "barcode")
ctrl_10_merged_metadata <- column_to_rownames(ctrl_10_merged_metadata, var = "barcode")
ctrl_10@meta.data <- ctrl_10_merged_metadata

ctrl_10 <- subset(ctrl_10, subset = doublet == "FALSE")
ctrl_10@meta.data$cells <- "LV_CAV"
ctrl_10@meta.data$sample_ID <- "6460_KA_10"

#------------------------------------------------------------------------------
# Create the merged file with all samples
merged.sets <- merge(ctrl_1, y = c(ctrl_2, ctrl_3, ctrl_4, ctrl_5, ctrl_6, ctrl_7, ctrl_8, ctrl_9, ctrl_10), add.cell.ids = c("RV_control", "RV_control", "RV_CAV", "LV_CAV", "RV_control", "RV_CAV", "RV_CAV", "LV_CAV", "RV_CAV", "LV_CAV"), project = "CAV_project")

features <- SplitObject(merged.sets, split.by = "sample_ID")
features <- lapply(X = features,
                   FUN = SCTransform,
                   method = "glmGamPoi",
                   vars.to.regress = "percent.mt",
                   return.only.var.genes = FALSE)
var.features <- SelectIntegrationFeatures(object.list = features, nfeatures = 3000)
merged.sets.sct <- merge(x = features[[1]], y = features[2:length(features)], merge.data=TRUE)
VariableFeatures(merged.sets.sct) <- var.features

# Continuing Seurat processing
merged.sets.sct <- RunPCA(merged.sets.sct, verbose = TRUE)

# Batch effect correction using Harmony
merged.sets.sct <- RunHarmony(merged.sets.sct, assay.use="SCT", group.by.vars = "cells")

# Continuing Seurat processing
merged.sets.sct <- RunUMAP(merged.sets.sct, reduction = "harmony", dims = 1:30)
merged.sets.sct <- FindNeighbors(merged.sets.sct, reduction = "harmony", dims = 1:30) %>% FindClusters(resolution = 0.4)

# Visualization of clusters
DimPlot(merged.sets.sct, group.by = "seurat_clusters", label = TRUE)

# View by sample label
DimPlot(merged.sets.sct, split.by = "cells", label = TRUE, ncol = 3)

# View by sample ID
DimPlot(merged.sets.sct, split.by = "sample_ID", label = TRUE, ncol = 5)

# Identification of cell types
# Cardiomyocytes
FeaturePlot(merged.sets.sct, features = c("TTN", "MLIP", "SORBS2", "FHL2"))

# Leukocytes
FeaturePlot(merged.sets.sct, features = c("PTPRC"))

# Macrophages
FeaturePlot(merged.sets.sct, features = c("CD163", "MRC1", "RBP1", "F13A1"))

# Mast cells
FeaturePlot(merged.sets.sct, features = c("MS4A2", "HPGDS", "HDC"))

# Lymphocyte
FeaturePlot(merged.sets.sct, features = c("SKAP1", "CD53", "PARP8", "PTPRC"))

# B cells
FeaturePlot(merged.sets.sct, features = c("MS4A1"))

# T cells
FeaturePlot(merged.sets.sct, features = c("CD4", "CD8A"))

# NK cells
FeaturePlot(merged.sets.sct, features = c("NCAM1", "CD160", "GNLY", "NKG7"))

# Fibroblasts
FeaturePlot(merged.sets.sct, features = c("DCN", "NEGR1", "ABCA6", "LAMA2"))

# Pericyte
FeaturePlot(merged.sets.sct, features = c("PDGFRB", "RGS5", "MCAM", "ACTA2"))

# Smooth muscle
FeaturePlot(merged.sets.sct, features = c("MYH11", "TAGLN", "ACTA2"))

# Neuronal
FeaturePlot(merged.sets.sct, features = c("NRXN1", "NCAM2"))

# Adipocytes
FeaturePlot(merged.sets.sct, features = c("ADIPOQ"))

# Endothelial cells
FeaturePlot(merged.sets.sct, features = c("VWF", "PECAM1", "LDB2", "ANO2"))

# Endocardial
FeaturePlot(merged.sets.sct, features = c("VWF", "PECAM1", "BMX", "NPR3"))

# Plasma cells
FeaturePlot(merged.sets.sct, features = c("SDC1", "IGKC"))

# Initial cluster labeling
updated.cluster.ids <- c("endothelial cells", "cardiomyocytes", "monocytes",
                         "CD8 T cells", "fibroblasts", "pericytes", "CD4 T cells",
                         "cardiomyocytes", "Tregs", "cardiomyocytes",
                         "NK cells", "endothelial cells", "B cells", "endocardial cells",
                         "VSMCs", "cardiomyocytes", "CD4 T cells", "endothelial cells",
                         "endothelial cells", "neuronal", "plasma cells", "dendritic cells", "cardiomyocytes", "mast cells", "endothelial cells")

names(updated.cluster.ids) <- levels(merged.sets.sct)
merged.sets.sct <- RenameIdents(merged.sets.sct, updated.cluster.ids)
DimPlot(merged.sets.sct, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend() + ggtitle("Annotated Clusters - Feb 27, 2022")

# Addressing monocyte/macrophage population using scGate data
library(scGate)

# Identifying the macrophage subpopulation
# Immune cells
immune <- subset(merged.sets.sct, idents = c("plasma cell", "NK cell", "macrophage", "mast cell", "B cell", "leukocyte"))
DimPlot(immune, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend() + ggtitle("Subset of Immune Cells")

immune@meta.data$cell_type <-  Idents(immune)
DefaultAssay(immune) <- "RNA"

# scGate analysis for macrophages
macrophage_scGate_model  <- gating_model(name = "macrophage", signature = c("CD163+", "MRC1+"))

immune <- scGate(data = immune, model = macrophage_scGate_model, assay = DefaultAssay(immune))
DimPlot(immune, group.by = "is.pure")
table(immune$is.pure)

macro_scGate <- subset(immune, subset = is.pure == "Pure")

# Pulling the macrophage cluster data into Seurat object
macro <- macro_scGate
Idents(macro) <- "macrophage"
macro_meta <- macro@meta.data
macro_meta$clusters <- "macrophage"
macro_meta <- dplyr::select(macro_meta, c(-is.pure, -is.pure.level1, -macrophage_UCell))

merged.sets.sct@meta.data$clusters <- Idents(merged.sets.sct)
merged_meta.data <- merged.sets.sct@meta.data
merged_meta.data <- rownames_to_column(merged_meta.data, var = "barcode")
macro_meta <- rownames_to_column(macro_meta, var = "barcode")
macro_meta <- dplyr::select(macro_meta, barcode, clusters)

test <- merged_meta.data %>% 
  mutate(coalesce(macro_meta$clusters[match(merged_meta.data$barcode, macro_meta$barcode)], merged_meta.data$clusters))
class(test$`coalesce(...)`)
class(test$clusters)
as.factor(test$`coalesce(...)`)

test$clusters <- test$`coalesce(...)`
test <- select(test, barcode:clusters)
test <- column_to_rownames(test, var = "barcode")
merged.sets.sct@meta.data <- test

Idents(merged.sets.sct) <- merged.sets.sct@meta.data$clusters

# Visualizing final clustered Seurat object including macrophages
DimPlot(merged.sets.sct, reduction = "umap", label = TRUE)
