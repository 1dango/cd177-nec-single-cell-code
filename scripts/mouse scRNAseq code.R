############################################################
## scRNA-seq analysis of mouse intestinal neutrophils
## R v4.2 | Seurat v4.3
############################################################

library(Seurat)
library(dplyr)
library(harmony)
library(ggplot2)
library(ggpubr)

set.seed(1234)

############################
# 1. Load Cell Ranger output
############################

NEC <- Read10X(data.dir = "NEC/filtered_feature_bc_matrix/")
CTRL <- Read10X(data.dir = "CTRL/filtered_feature_bc_matrix/")

NEC <- CreateSeuratObject(NEC, min.cells = 3, min.features = 200)
CTRL <- CreateSeuratObject(CTRL, min.cells = 3, min.features = 200)

NEC$group <- "NEC"
CTRL$group <- "CTRL"
NEC$sample <- "NEC"
CTRL$sample <- "CTRL"

############################
# 2. Merge datasets
############################

obj <- merge(NEC, CTRL, add.cell.ids = c("NEC", "CTRL"))

############################
# 3. Quality control
############################

obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = "^Rp[sl]")

obj <- subset(
  obj,
  subset =
    nFeature_RNA > 200 &
    nFeature_RNA < 8000 &
    nCount_RNA > 300 &
    nCount_RNA < 60000 &
    percent.mt < 25
)

############################
# 4. Normalization and scaling
############################

obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, nfeatures = 2000)
obj <- ScaleData(obj, vars.to.regress = c("percent.mt", "percent.ribo"))

############################
# 5. Dimensional reduction & integration
############################

obj <- RunPCA(obj)
obj <- RunHarmony(obj, group.by.vars = "sample")

obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:10)
obj <- FindClusters(obj, resolution = 0.2)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:10)

############################
# 6. Cell type annotation
############################

new.cluster.ids <- c(
  "enterocytes", "stem_cells", "enterocytes",
  "neutrophils", "stem_cells", "fibroblasts",
  "macrophages", "goblet_cells", "enterocytes"
)

names(new.cluster.ids) <- levels(obj)
obj <- RenameIdents(obj, new.cluster.ids)
obj$celltype <- Idents(obj)

############################
# 7. Neutrophil subclustering
############################

neu <- subset(obj, idents = "neutrophils")

neu <- FindNeighbors(neu, reduction = "harmony", dims = 1:10)
neu <- FindClusters(neu, resolution = 0.1)
neu <- RunUMAP(neu, reduction = "harmony", dims = 1:10)

neu <- RenameIdents(neu, c("Il11b_neu", "Cd177_neu"))
neu$celltype <- Idents(neu)

############################
# 8. Differential expression
############################

markers <- FindAllMarkers(
  neu,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

############################
# 9. Gene expression comparison (violin plots)
############################

genes.use <- c("Tpm4", "Actn1", "Wdr1", "Myh9")

VlnPlot(
  neu,
  features = genes.use,
  group.by = "celltype",
  pt.size = 0.1
)

############################
# 10. Statistical testing
############################

pvals <- sapply(genes.use, function(g) {
  x <- FetchData(neu, g)[neu$celltype == "Il11b_neu", 1]
  y <- FetchData(neu, g)[neu$celltype == "Cd177_neu", 1]
  t.test(x, y)$p.value
})

