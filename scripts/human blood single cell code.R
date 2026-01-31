#---blood----
# run_all_blood.R  (English comments version)
rm(list = ls())
setwd("/home/dell/Project/General/NEC/single_cell/blood")   # change if necessary
raw_data_dir <- "raw"
dir.create("results", showWarnings = FALSE)

## Load required libraries
library(Seurat)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(DoubletFinder)
library(Harmony)
library(clusterProfiler)
library(org.Hs.eg.db)
library(pheatmap)
library(RColorBrewer)
library(cowplot)

theme_nec <- theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = "bold"))

## 1. Load data
samples <- list.dirs(raw_data_dir, full.names = FALSE, recursive = FALSE)
obj.list <- map(samples, function(s) {
  counts <- Read10X(file.path(raw_data_dir, s))
  CreateSeuratObject(counts, project = s, min.features = 200, min.cells = 3)
})

## 2. Merge objects
combined <- merge(obj.list[[1]], y = obj.list[-1], add.cell.ids = samples)

## 3. Quality control
eryth.genes <- c("HBA1", "HBA2", "HBB", "HBG1", "HBG2", "HBZ", "HBM", "AHSP")
combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^MT-")
combined[["percent.ery"]] <- PercentageFeatureSet(combined, features = eryth.genes)
combined <- subset(combined, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 &
                     percent.mt < 5 & percent.ery < 10)

## 4. Doublet removal
sweep.res <- paramSweep(combined, num.cores = 1)
bcmvn <- find.pK(sweep.res)
pk.opt <- as.numeric(as.matrix(bcmvn)[which.max(bcmvn$BCmetric), "pK"])
combined <- doubletFinder(combined, pK = pk.opt, expected.doublets = 0.05)
combined <- subset(combined, subset = DF.classifications == "Singlet")

## 5. Harmony integration
combined <- NormalizeData(combined, normalization.method = "LogNormalize", scale.factor = 10000)
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)
combined <- ScaleData(combined)
combined <- RunPCA(combined, features = VariableFeatures(combined))
combined <- RunHarmony(combined, group.by.vars = "orig.ident",
                       assay.use = "RNA", theta = 2, lambda = 1, max.iter.harmony = 20)
combined <- FindNeighbors(combined, reduction = "harmony", dims = 1:30)
combined <- FindClusters(combined, resolution = 0.5)
combined <- RunUMAP(combined, reduction = "harmony", dims = 1:30)

## 6. Add group metadata
group_vector <- case_when(
  all_cells@meta.data$orig.ident %in% c( “NEC-1”,"NEC-2", "NEC-3") ~ "NEC",
  all_cells@meta.data$orig.ident %in% c("Ctrl-1", "Ctrl-2") ~ "CTRL",
  TRUE ~ as.character(all_cells@meta.data$orig.ident)  # 确保返回字符型
)
all_cells <- AddMetaData(all_cells, metadata = group_vector, col.name = "group")
## 7. Rename cell clusters (blood annotation)
new.ids <- c(
  'T cells', 'Neutrophil', 'Monocyte', 'Monocyte', 'Macrophages', 'T cells', "T cells",
  'B cells', 'Natural killer cell', 'Neutrophil', 'Macrophages', "T cells", 'Neutrophil', "T cells", "Myeloid stem cells",
  "T cells", "T cells", "Macrophages", "Platelet", "Natural killer cell", "T cells", "Dendritic cells",
  "Hematopoietic stem cell", "B cells", "Neutrophil"
)
names(new.ids) <- levels(combined)
combined <- RenameIdents(combined, new.ids)

## 8. KEGG pathway scoring
kegg.genes <- read_csv("hsahsa04610&hsa04611.csv", col_names = "gene") %>% pull(gene)
combined <- AddModuleScore(combined, features = list(kegg.genes),
                           ctrl = 100, name = "KEGG_coag_platelet")

## 9. Violin plot of pathway score
p.vln <- VlnPlot(combined, features = "KEGG_coag_platelet1", split.by = "group",
                 group.by = "ident", pt.size = 0, split.plot = TRUE) + theme_nec
ggsave("results/09_vln_kegg_coagPlatelet_blood.pdf", p.vln, width = 8, height = 5)

## 10. Heatmap with gene_order
gene_order <- c("VASP", "MYL12A", "TBXAS1", "SYK", "FERMT3", "CD55",
                "ITGAX", "ITGAM", "ITGB2", "CR1", "PLAUR", "PIK3R5",
                "RHOA", "ROCK1", "MYL12B", "GNAQ", "APBB1IP", "TLN1",
                "GNAI2", "ACTG1", "ACTB", "FCGR2A", "FCER1G", "LYN")
avg.exp <- AverageExpression(combined, group.by = "group", assays = "RNA")$RNA %>%
  select(all_of(gene_order)) %>% as.matrix()
tissue_T <- t(scale(t(avg.exp)))
color_palette <- colorRampPalette(c("#4E76B4", "#EEF0EF", "#EA4C5A"))(100)
pheatmap(tissue_T,
         scale = "column", cluster_rows = FALSE, cluster_cols = FALSE,
         color = color_palette, show_rownames = TRUE, show_colnames = TRUE,
         cellwidth = 10, cellheight = 10, fontsize_row = 8, fontsize_col = 8,
         filename = "results/10_heatmap_geneOrder_blood.pdf")

## 11. Re-cluster neutrophils
neu <- subset(combined, idents = "Neutrophil")
neu <- NormalizeData(neu, normalization.method = "LogNormalize", scale.factor = 10000)
neu <- FindVariableFeatures(neu, selection.method = "vst", nfeatures = 2000)
neu <- ScaleData(neu)
neu <- RunPCA(neu, features = VariableFeatures(neu))
neu <- RunHarmony(neu, group.by.vars = "orig.ident", theta = 2, lambda = 1, max.iter.harmony = 20)
neu <- FindNeighbors(neu, reduction = "harmony", dims = 1:20)
neu <- FindClusters(neu, resolution = 0.5)
neu <- RunUMAP(neu, reduction = "harmony", dims = 1:20)

## 12. Feature plots
# CD177 across all neutrophils
p.cd177 <- FeaturePlot(neu, features = "CD177", split.by = "group",
                       pt.size = 0.5, cols = c("lightgrey", "red")) +
  theme_nec + ggtitle("CD177 expression in neutrophils (blood)")
ggsave("results/11a_feature_cd177_blood.pdf", p.cd177, width = 6, height = 4)

# NEC vs CTRL for selected genes
target.genes <- c("CD177", "IL1B", "FOS")
p.split <- FeaturePlot(neu, features = target.genes, split.by = "group",
                       pt.size = 0.5, cols = c("lightgrey", "blue")) +
  theme_nec + plot_layout(guides = "collect")
ggsave("results/11b_feature_splitNEC_CTRL_blood.pdf", p.split, width = 7, height = 4)

## 13. KEGG enrichment (NEC vs CTRL DEGs, Interferon-hi background)
degs <- FindMarkers(neu, ident.1 = "NEC", ident.2 = "CTRL",
                    min.pct = 0.1, logfc.threshold = 0)
write_csv(degs %>% rownames_to_column("gene"), "results/12_nec_vs_ctrl_degs_blood.csv")

fdegs <- read_csv("Interferon_hi_Neu.csv")          # must contain SYMBOL column
need_fDEG <- fdegs[, c(1, 2, 3)]
colnames(need_fDEG) <- c("SYMBOL", "pvalue", "log2FoldChange")
df <- bitr(need_fDEG$SYMBOL, fromType = "SYMBOL",
           toType = "ENTREZID", OrgDb = org.Hs.eg.db)
need_fDEG <- merge(need_fDEG, df, by = "SYMBOL")

geneList <- need_fDEG$log2FoldChange
names(geneList) <- need_fDEG$ENTREZID
geneList <- sort(geneList, decreasing = TRUE)

kk <- enrichKEGG(gene = geneList,
                 organism = "hsa",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.2,
                 universe = names(geneList))   # background = Interferon-hi genes

if (!is.null(kk) && nrow(kk@result) > 0) {
  kk.df <- as.data.frame(kk) %>% arrange(pvalue)
  p.kegg <- ggplot(kk.df, aes(x = Count, y = reorder(Description, -pvalue),
                              fill = -log10(pvalue))) +
    geom_bar(stat = "identity", width = 0.7, alpha = 0.8) +
    scale_fill_distiller(palette = "YlOrRd", direction = 1) +
    labs(x = "Gene count", y = NULL,
         title = "KEGG enrichment (NEC vs CTRL, Interferon-hi background)") +
    theme_nec
  ggsave("results/13_kegg_nec_vs_ctrl_interferon_blood.pdf", p.kegg, width = 6, height = 4)
}

## 14. Generate neu_dot.csv & dotplot
dot.ctrl <- FindMarkers(neu, ident.1 = "CTRL", group.by = "group",
                        min.pct = 0.1, logfc.threshold = 0) %>%
  rownames_to_column("gene") %>% filter(avg_log2FC > 0) %>%
  slice_head(n = 50) %>% mutate(cluster = "CTRL")

dot.nec <- FindMarkers(neu, ident.1 = "NEC", group.by = "group",
                       min.pct = 0.1, logfc.threshold = 0) %>%
  rownames_to_column("gene") %>% filter(avg_log2FC > 0) %>%
  slice_head(n = 50) %>% mutate(cluster = "NEC")

dot.df <- bind_rows(dot.ctrl, dot.nec) %>%
  mutate(pct = map2_dbl(gene, cluster, function(g, cl) {
    100 * mean(neu@assays$RNA@data[g, neu$group == cl] > 0)
  })) %>%
  select(gene, cluster, avg_log2FC = avg_log2FC, pct)

write_csv(dot.df, "results/14_neu_dot_blood.csv")

dot.df$gene <- factor(dot.df$gene, levels = unique(dot.df$gene))
dot.df$cluster <- factor(dot.df$cluster, levels = c("CTRL", "NEC"))
p.dot <- ggplot(dot.df, aes(x = gene, y = cluster)) +
  geom_point(aes(size = pct, fill = avg_log2FC), shape = 21,
             color = "black", stroke = 0.5) +
  scale_size_area(name = "Percent\nexpressed", max_size = 8, limits = c(0, 100)) +
  scale_fill_gradient2(low = "#4E76B4", mid = "#EEF0EF", high = "#EA4C5A",
                       limits = c(-1, 2), name = "Avg log2FC") +
  theme_nec +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
ggsave("results/15_neu_dotplot_blood.pdf", p.dot, width = 7, height = 4)

## 15. Save final objects
saveRDS(combined, "results/16_final_seurat_blood.rds")
saveRDS(neu, "results/16_final_neu_blood.rds")