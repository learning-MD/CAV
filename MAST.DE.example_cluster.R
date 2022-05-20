library(Seurat)
library(tidyverse)
library(MAST)

# Endothelial cells
endo <- subset(merged.sets.sct, idents = c("endothelial cells"))

# Converting Seurat object into SingleCellExperiment
sce <- as.SingleCellExperiment(DietSeurat(endo, graphs = c("pca", "umap")))

# Log2 transform the raw data
assay(sce, "logcounts") <- log2(counts(sce) + 1)

# Convert the SingleCellExperiment object to SingleCellAssay required by MAST
sca <- MAST::SceToSingleCellAssay(sce, class = "SingleCellAssay")

# Scale gene detection rate
colData(sca)$nFeature_RNA = scale(colData(sca)$nFeature_RNA)

cdr2 <- colSums(assay(sca)>0)
colData(sca)$cngeneson <- scale(cdr2)

cond <- factor(colData(sca)$cells) #Variable of interest
cond <- relevel(cond, "RV_control") #control
colData(sca)$condition <- cond

# Filtering for genes with non-zero expression in >5% of nuclei
expressed_genes <-MAST::freq(sca) >0.05
sca <- sca[expressed_genes,]

# Running MAST
zlmCond <- suppressMessages(MAST::zlm(~condition + cngeneson + (1 | sample_ID), exprs_value = "logcounts", sca, method='glmer',ebayes = FALSE,strictConvergence = FALSE, parallel=TRUE))

summaryCond <- summary(zlmCond, doLRT='conditionRV_CAV')

summaryDt <- summaryCond$datatable

fcHurdle <- merge(summaryDt[contrast=='conditionRV_CAV' & component=='H',.(primerid, `Pr(>Chisq)`)],
                  summaryDt[contrast=='conditionRV_CAV' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid')

fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
View(fcHurdle)

write_csv(fcHurdle, "endothelial_cells_MAST.csv")
