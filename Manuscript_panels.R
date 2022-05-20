# Loading packages
library(Seurat)
library(tidyverse)
library(EnhancedVolcano)

# Loading data
merged.sets.sct <- readRDS()
chimerism <- readRDS()

# Figure 1B
DimPlot(merged.sets.sct)

# Figure 1E
endo_deg <- read_csv("endothelial_cells_MAST.csv")
p1 <- EnhancedVolcano(endo_deg,
                      lab = endo_deg$primerid,
                      x = "coef",
                      y = "fdr",
                      title = "Endothelial Cells - RV CAV vs RV Control",
                      subtitle = NULL,
                      ylab = "-Log10 FDR",
                      FCcutoff = 0.1,
                      pCutoff = 0.05,
                      col = c("gray", "gray", "gray", "red"),
                      labSize = 3,
                      legendPosition = "none",
                      xlim = c(-1.25, 1.25),
                      ylim = c(0, 4))
p1

# Figure 1G
macro_deg <- read_csv("macrophages_MAST.csv")
p2 <- EnhancedVolcano(macro_deg,
                      lab = macro_deg$primerid,
                      x = "coef",
                      y = "fdr",
                      title = "Macrophages - RV CAV vs RV Control",
                      subtitle = NULL,
                      ylab = "-Log10 FDR",
                      FCcutoff = 0.1,
                      pCutoff = 0.05,
                      col = c("gray", "gray", "gray", "red"),
                      labSize = 3,
                      legendPosition = "none",
                      xlim = c(-1.25, 1.25),
                      ylim = c(0, 4))
p2

# Figure 1H
DimPlot(chimerism, group.by = "origin")

# Figure 1I
library(scGate)

endo <- subset(chimerism, idents = "endothelial cells")
DimPlot(endo, split.by = "origin")
DefaultAssay(endo) <- "RNA"

my_scGate_model <- gating_model(name = "EndoMT", signature = c("SERPINE1", "VIM", "FN1"))
seurat_object <- scGate(data = endo, model = my_scGate_model)
DimPlot(seurat_object, group.by = "is.pure", split.by = "origin") +
  ggtitle("EndoMT (SERPINE1+, VIM+, FN1+)") +
  theme(plot.title = element_text(size = 12, face = "bold"))

# Figure 1J
library(Nebulosa)
macro <- subset(chimerism, idents = "monocytes")
plot_density(macro, features = c("CCR2", "MRC1"))
