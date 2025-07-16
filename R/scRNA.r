setwd
library(Seurat)
library(harmony)
library(tidyverse)
library(sctransform)
library(cols4all)
library(clustree)
library(ggpubr)
library(stringr)
library(tidydr)
library(openxlsx)
library(data.table)
library(reshape2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(clusterProfiler)
library(pheatmap)
library(ComplexHeatmap)
library(GSVA)
library(GSEABase)
library(fgsea)
library(corrplot)
library(colorspace)
library(survival)
library(survminer)
library(maftools)
library(vegan)
library(forcats)
library(ggpubr)
library(ggsci)
library(ggplot2)
library(rstatix)
library(ggcor)
library(ggstance)
library(ggplot2)
library(magrittr)
library(gtools)
library(stringr)
library(tidyverse)
library(patchwork)
library(data.table)
library(RColorBrewer)
library(cowplot)
library(ggpubr)
dir <- dir("/home/pub252/users/wangxiaonan/20250414_OV/00.rawdata/GSE184880/")
samples_name = dir
scRNAlist <- list()
for(i in 1:length(dir)){
  counts <- Read10X(data.dir=dir[i])
  scRNAlist[[i]] <- CreateSeuratObject(counts, project=samples_name[i], min.cells=3, min.features=200)
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]], add.cell.id=samples_name[i])   
  if(T){    
    scRNAlist[[i]][["percent.mito"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern="^MT-") 
  }
}
names(scRNAlist) <- samples_name
scRNA <- merge(scRNAlist[[1]], scRNAlist[2:length(scRNAlist)])

scRNA = subset(scRNA, subset=nFeature_RNA>200&nFeature_RNA<6000&percent.mito<10)  # 过滤细胞 16526
colors = c("orange", "chartreuse", "cyan","red" ,"deepskyblue", "brown",  "darkviolet",  "hotpink","magenta", "mediumblue","mediumturquoise")
p = VlnPlot(scRNA, features=c("nFeature_RNA"), pt.size=0, cols=colors)
ggsave("nFeature_RNA_violin.pdf", p, width=6, height=6)

scRNA = SCTransform(scRNA, verbose=FALSE)
scRNA = RunPCA(scRNA, verbose=FALSE)
colnames(scRNA@meta.data)[1] = "Sample"
scRNA = RunHarmony(scRNA, group.by.vars="Sample", max.iter.harmony=50, lambda=0.5, assay.use="SCT")
ElbowPlot(scRNA)
scRNA <- FindNeighbors(scRNA, dims=1:50, reduction="harmony")
scRNA <- RunUMAP(scRNA, dims=1:50, reduction="harmony")
###############################################################
custome_theme_1=theme(axis.text.y=element_text(family="Times",face="plain")
                      ,axis.text.x=element_text(family="Times",face="plain")
                      ,plot.title = element_text(hjust = 0.5,family="Times",face="plain")
                      ,axis.title.x=element_text(family="Times",face="plain")
                      ,axis.title.y=element_text(family="Times",face="plain")
                      ,legend.title = element_text(family="Times",face="plain")
                      ,legend.text = element_text(family="Times",face="plain"))
sc_umap = DimPlot(scRNA,cols=colors,group.by='Sample',
                  reduction="umap",
                  label = F, 
                  pt.size = 0.2,
                  label.size = 5) +
  theme(axis.line = element_line(size=0.1, colour = "black"), 
        #axis.text = element_blank(), 
        #axis.title = element_blank(),
        axis.ticks = element_blank()
  ) +ggtitle('')
sc_umap
# LabelClusters(sc_umap,id = 'cell_type',family = 'Times',size = 6,fontface = 'bold',color = 'red')
sc_umap_Sample = DimPlot(scRNA,cols=colors,group.by='Sample',
                         reduction="umap",
                         label = "F", 
                         pt.size = 0.2,
                         label.size = 0) +
  theme(axis.line = element_line(size=0.1, colour = "black"), 
        #axis.text = element_blank(), 
        #axis.title = element_blank(),
        axis.ticks = element_blank()
  ) +ggtitle('')

sc_umap_Sample
p1=sc_umap+custome_theme_1
p2=sc_umap_Sample+custome_theme_1
p1=ggarrange(AugmentPlot(p1,dpi = 300), ncol = 1, nrow = 1,
             legend ='right',legend.grob=get_legend(p1),
             font.label = list(size = 12, face = "bold",family ='Times')) # 设置标签字体样式
ggsave(filename = 'sample_umap.pdf',plot = p1,he=8,width = 9)
#########################################################################
scRNA <- FindClusters(scRNA, resolution = seq(0.05,0.5,by=0.05))
clustree(scRNA)
mydata <- FindClusters(scRNA, resolution=0.05)
UMAPPlot(mydata, pt.size=1, label=T, cols=colors)+NoLegend()
markers <- FindAllMarkers(mydata, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.30)
write.table(markers, "subcluster_DEG.txt", col.names=T, row.names=T, quote=F, sep="\t")

cell_label = c("T cell","Epithelial cell","Fibroblast","Macrophage","Plasma cell","MKI67+ progenitor cell","B cell","Plasmacytoid dendritic cell")
#C0:T cell:CCL5	NKG7	GNLY	CD2
#C1:Epithelial cell:WFDC2	EPCAM	KRT18	KRT8
#C2:Fibroblast:LUM	DCN	VCAN
#C3:Macrophage:LYZ	C1QA	C1QB	C1QC
#C4:Plasma cell:IGHG1	IGHG3	IGHG4	CD79A
#C5:MKI67+ progenitor cell:MKI67 TOP2A
#C6:B cell:MS4A1	CD79A	CD22
#C7:Plasmacytoid dendritic cell:IL3RA	LILRA4	PLD4
names(cell_label) <- levels(mydata)
mydata <- RenameIdents(mydata, cell_label)
mydata[["cell_type"]] = Idents(mydata)
mydata = subset(mydata, cell_type %in% c("T cell","Epithelial cell","Fibroblast","Macrophage","Plasma cell","MKI67+ progenitor cell","B cell","Plasmacytoid dendritic cell"))
###############################################################
custome_theme_1=theme(axis.text.y=element_text(family="Times",face="plain")
                      ,axis.text.x=element_text(family="Times",face="plain")
                      ,plot.title = element_text(hjust = 0.5,family="Times",face="plain")
                      ,axis.title.x=element_text(family="Times",face="plain")
                      ,axis.title.y=element_text(family="Times",face="plain")
                      ,legend.title = element_text(family="Times",face="plain")
                      ,legend.text = element_text(family="Times",face="plain"))
sc_umap = DimPlot(mydata,cols=colors,group.by='cell_type',
                  reduction="umap",
                  label = F, 
                  pt.size = 0.2,
                  label.size = 5) +
  theme(axis.line = element_line(size=0.1, colour = "black"), 
        #axis.text = element_blank(), 
        #axis.title = element_blank(),
        axis.ticks = element_blank()
  ) +ggtitle('')
sc_umap
# LabelClusters(sc_umap,id = 'cell_type',family = 'Times',size = 6,fontface = 'bold',color = 'red')
sc_umap_Sample = DimPlot(mydata,cols=colors,group.by='cell_type',
                         reduction="umap",
                         label = "F", 
                         pt.size = 0.2,
                         label.size = 0) +
  theme(axis.line = element_line(size=0.1, colour = "black"), 
        #axis.text = element_blank(), 
        #axis.title = element_blank(),
        axis.ticks = element_blank()
  ) +ggtitle('')

sc_umap_Sample
p1=sc_umap+custome_theme_1
p2=sc_umap_Sample+custome_theme_1
p1=ggarrange(AugmentPlot(p1,dpi = 300), ncol = 1, nrow = 1,
             legend ='right',legend.grob=get_legend(p1),
             font.label = list(size = 12, face = "bold",family ='Times'))
ggsave(filename = 'cell_umap.pdf',plot = p1,he=8,width = 9)
#########################################################################
saveRDS(mydata, "scRNA_raw.rds")

genes = c("CCL5","NKG7","GNLY","CD2","WFDC2","EPCAM","KRT18","KRT8","LUM","DCN","VCAN","LYZ","C1QA","C1QB","C1QC",
          "IGHG1","IGHG4","IGHG3","CD79A","MKI67","TOP2A","MS4A1","CD22","IL3RA","LILRA4","PLD4")
p = DotPlot(mydata, features=genes)+coord_flip()+scale_color_distiller(palette="RdYlBu")+
  theme(text=element_text(family="Times"), axis.text.x=element_text(angle=30, hjust=1, size=12, face="bold"), axis.text.y=element_text(face="bold", size=12), axis.title.x=element_blank(), axis.title.y=element_blank())
ggsave("marker_dotplot.pdf", p, width=7, height=6)

p = DoHeatmap(mydata, features=genes, group.colors=colors, label=FALSE)+scale_fill_gradientn(colors=c("white", "snow", "firebrick3"))+theme(text=element_text(family="Times"))
ggsave("marker_heatmap.pdf", p, width=10, height=5)

bar = mydata@meta.data %>% group_by(Sample, cell_type) %>% count()
bar$cell_type = factor(bar$cell_type, levels=cell_label)
p = ggplot(data=bar, aes(x=Sample, y=n, fill=cell_type))+ 
  geom_bar(stat="identity", position=position_fill())+
  scale_fill_manual(values=colors)+theme_bw()+
  theme(text=element_text(family="Times"), axis.title.x=element_blank(), axis.title.y=element_blank(), 
        axis.text.x=element_text(face="bold", size=8, angle=15, hjust=1), axis.text.y=element_text(face="bold", size=12),
        legend.text=element_text(face="bold", size=8), legend.title=element_blank(), legend.position="right")
ggsave("barplot_cellType.pdf", p, width=8, height=6)

cell_DEGs.list=split(x=markers,f=markers$cluster)
cell_DEGs.list=sapply(cell_DEGs.list, function(x){subset(x,select='gene',drop=TRUE)})
library(topGO)
library(clusterProfiler)
library(Rgraphviz)
library(org.Hs.eg.db)
library(ggplot2)
library(cowplot)
library(ggpubr)
subcluster.name=names(cell_DEGs.list)
subcluster.name=subcluster.name[order(subcluster.name)]
subcluster.name=c("2")
p=list()
for (s in subcluster.name) {
  erich.go.BP = enrichGO(gene =  cell_DEGs.list[[s]],OrgDb = org.Hs.eg.db,keyType = "SYMBOL",ont = "BP",pvalueCutoff = 0.05)
  erich.go.BP.res=erich.go.BP@result
  write.csv(erich.go.BP.res,paste0('./',s,'_enrichment.csv'))
  erich.go.BP.res=erich.go.BP.res[erich.go.BP.res$p.adjust<0.05,c('Description','GeneRatio','Count','p.adjust')]
  erich.go.BP.res=erich.go.BP.res %>% slice_min(n =10, order_by = p.adjust)
  p[[which(subcluster.name==s)]]=ggplot(data = erich.go.BP.res,
                                        mapping = aes(x=Count,y=reorder(Description,Count), fill = -log10(p.adjust))) +
    geom_bar(stat="identity")+ theme_bw()+
    scale_fill_gradient(low="yellow",high =c("red")[which(subcluster.name==s)])+
    labs(y=NULL,title = "subcluster ")+ggtitle(s)+
    theme(text = element_text(family = 'Times',size = 14))+
    scale_y_discrete(labels=function(y)str_wrap(y,width = 35))
  
}
length(p)
subcluster.enrichment.fig <- plot_grid(plotlist = p, ncol = 3, nrow = 2)
pdf('subcluster_enrichment.pdf',height = 20,width = 20,onefile = F)
subcluster.enrichment.fig
dev.off()


