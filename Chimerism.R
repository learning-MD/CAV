library(Seurat)
library(tidyverse)
library(patchwork)
library(DropletUtils)
library(glmGamPoi)
library(sctransform)
library(harmony)
library(SeuratDisk)
memory.limit(10000000000)

# To load sample
merged.sets.sct <- readRDS("merged.sets.sct.rds")

# For chimerism analysis, use 4, 7, 8, 9, 10
chimerism <- subset(merged.sets.sct,
                    subset = sample_ID == c("6460_KA_4", "6460_KA_7", "6460_KA_8", "6460_KA_9", "6460_KA_10"))

# Subsetting individual samples
KA_4 <- subset(chimerism, subset = sample_ID == "6460_KA_4")
KA_4.meta.data <- KA_4@meta.data
KA_4.meta.data <- KA_4.meta.data %>%
  mutate(origin = str_replace_all(KA_4.meta.data$assignment,
                                  c("1" = "donor", "0" = "recipient")))
KA_4.meta.data <- rownames_to_column(KA_4.meta.data, var = "barcode")
KA_4.meta.data <- dplyr::select(KA_4.meta.data, barcode, origin)

KA_7 <- subset(chimerism, subset = sample_ID == "6460_KA_7")
KA_7.meta.data <- KA_7@meta.data
KA_7.meta.data <- KA_7.meta.data %>%
  mutate(origin = str_replace_all(KA_7.meta.data$assignment,
                                  c("1" = "donor", "0" = "recipient")))
KA_7.meta.data <- rownames_to_column(KA_7.meta.data, var = "barcode")
KA_7.meta.data <- dplyr::select(KA_7.meta.data, barcode, origin)

KA_8 <- subset(chimerism, subset = sample_ID == "6460_KA_8")
KA_8.meta.data <- KA_8@meta.data
KA_8.meta.data <- KA_8.meta.data %>%
  mutate(origin = str_replace_all(KA_8.meta.data$assignment,
                                  c("0" = "donor", "1" = "recipient")))
KA_8.meta.data <- rownames_to_column(KA_8.meta.data, var = "barcode")
KA_8.meta.data <- dplyr::select(KA_8.meta.data, barcode, origin)

KA_9 <- subset(chimerism, subset = sample_ID == "6460_KA_9")
KA_9.meta.data <- KA_9@meta.data
KA_9.meta.data <- KA_9.meta.data %>%
  mutate(origin = str_replace_all(KA_9.meta.data$assignment,
                                  c("0" = "donor", "1" = "recipient")))
KA_9.meta.data <- rownames_to_column(KA_9.meta.data, var = "barcode")
KA_9.meta.data <- dplyr::select(KA_9.meta.data, barcode, origin)

KA_10 <- subset(chimerism, subset = sample_ID == "6460_KA_10")
KA_10.meta.data <- KA_10@meta.data
KA_10.meta.data <- KA_10.meta.data %>%
  mutate(origin = str_replace_all(KA_10.meta.data$assignment,
                                  c("1" = "donor", "0" = "recipient")))
KA_10.meta.data <- rownames_to_column(KA_10.meta.data, var = "barcode")
KA_10.meta.data <- dplyr::select(KA_10.meta.data, barcode, origin)

# Stack together the metadata files
total.metadata <- full_join(KA_4.meta.data, KA_7.meta.data, by = c("barcode", "origin"))
total.metadata <- full_join(total.metadata, KA_8.meta.data, by = c("barcode", "origin"))
total.metadata <- full_join(total.metadata, KA_9.meta.data, by = c("barcode", "origin"))
total.metadata <- full_join(total.metadata, KA_10.meta.data, by = c("barcode", "origin"))

# Create chimerism metadata file
chimerism.meta.data <- chimerism@meta.data
chimerism.meta.data <- rownames_to_column(chimerism.meta.data, var = "barcode")

# Join chimerism metadata with the individual metadata
chimerism.meta.data <- left_join(chimerism.meta.data, total.metadata, by = "barcode")
chimerism.meta.data <- column_to_rownames(chimerism.meta.data, var = "barcode")
chimerism@meta.data <- chimerism.meta.data

# Subset only the "donor" and "recipient" cells
chimerism <- subset(chimerism, subset = origin == c("donor", "recipient"))

DimPlot(chimerism, reduction = "umap")
DimPlot(chimerism, reduction = "umap", split.by = "origin")
DimPlot(chimerism, reduction = "umap", split.by = "sample_ID")