setwd
library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)
library(igraph)
theme_set(theme_cowplot())
set.seed(123)
seurat_obj <- readRDS("scRNA_raw.rds")
DimPlot(seurat_obj, group.by = "cell_type", label = T) + umap_theme()+ ggtitle("GSE184880") + NoLegend()

seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction",
  fraction = 0.05,
  wgcna_name = "GSE184880" # hdWGCNA
)

seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("cell_type", "Sample"), 
  reduction = "harmony", 
  k = 25, 
  max_shared = 20, 
  ident.group = "cell_type"  
)

seurat_obj <- NormalizeMetacells(seurat_obj)
seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = "Fibroblast", 
  group.by = "cell_type", 
  assay = "SCT",
  slot = "data"
)

seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = "signed" 
)

plot_list <- PlotSoftPowers(seurat_obj)

a<-wrap_plots(plot_list, ncol=2)
ggsave("Soft Threshold.pdf", a, width=10, height=5)

seurat_obj <- ConstructNetwork(
  seurat_obj, soft_power = 12, 
  setDatExpr = FALSE,
  tom_name = 'Fibroblast' # TOM

pdf("Dendrogram.pdf", width = 8, height = 6)
PlotDendrogram(seurat_obj, main='Fibroblast hdWGCNA Dendrogram')
dev.off()

TOM <- GetTOM(seurat_obj)
TOM[1:4, 1:4]

seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(seurat_obj))

seurat_obj <- ModuleEigengenes(
  seurat_obj,
  group.by.vars = "Sample"  
)

seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = 'cell_type', group_name = 'Fibroblast'  
)

seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = "Fibroblast-M"
)

pdf("module.pdf", width = 15, height = 8)
PlotKMEs(seurat_obj, ncol = 5, n_hubs = 10)
dev.off()

modules <- GetModules(seurat_obj)

hub_df <- GetHubGenes(seurat_obj = seurat_obj, n_hubs = 10)

saveRDS(seurat_obj, file = 'hdWGCNA_object.rds')

seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 25, # topN hub genes
  method = "Seurat" ，Seurat、UCell
)

plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features = 'hMEs', MEs、hMEs、scores、average
  order = TRUE # order so the points with highest hMEs are on top
)

#  stitch together with patchwork
pdf("module_score.pdf", width = 15, height = 8)
wrap_plots(plot_list, ncol = 6)
dev.off()

pdf("module_cor.pdf", width = 8, height = 8)
ModuleCorrelogram(seurat_obj)
dev.off()

MEs <- GetMEs(seurat_obj, harmonized=TRUE)
mods <- colnames(MEs); mods <- mods[mods != 'grey']

# add hMEs to Seurat meta-data:
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)

p <- DotPlot(seurat_obj, features = mods, group.by = "cell_type")

p <- p + 
  #     coord_flip() + 
  RotatedAxis() + 
  scale_color_gradient2(high = "red", mid = "grey95", low = "blue")  

p
ggsave("module_cell.pdf", p, width=8, height=8)



ModuleNetworkPlot(seurat_obj = seurat_obj, mods = "Fibroblast-M1")
ModuleNetworkPlot(seurat_obj = seurat_obj, mods = "Fibroblast-M5")
ModuleNetworkPlot(seurat_obj = seurat_obj, mods = "Fibroblast-M8")
ModuleNetworkPlot(seurat_obj = seurat_obj, mods = "Fibroblast-M11")
ModuleNetworkPlot(seurat_obj = seurat_obj, mods = "Fibroblast-M7")

pdf("module_network.pdf", width = 8, height = 8)
HubGeneNetworkPlot(
  seurat_obj,
  n_hubs = 10, hub gene
  n_other = 0, gene
  edge_prop = 0.75, 
  mods = "all"
)
dev.off()

