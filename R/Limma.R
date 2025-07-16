setwd
library(limma)
library(pheatmap)
library(ggplot2)
library(data.table)
library(dplyr)
data_GEO <- read.csv("GSE66957.csv",header = T,check.names = F,row.names = 1,sep = ",")
Control_num<-12
Case_num<-57
logFoldChange <- 1
P.Value <- 0.01 
adj.P.Val <- 0.01

group<-factor(c(rep("Control",Control_num),rep("Case",Case_num)),levels=c('Control','Case')) 
design<-model.matrix(~0+group)
colnames(design)<-c("Control","Case")
fit<-lmFit(data_GEO,design)
contrast.matrix <- makeContrasts(Case-Control,levels=design)
fit2<-contrasts.fit(fit,contrast.matrix)
fit2<-eBayes(fit2)
options(digits = 4)
result_mRNA<-topTable(fit2,coef=1,adjust="BH",number=dim(data_GEO)[1])
write.table(result_mRNA,file="DEG_information.txt",sep="\t",quote=F)
diffSig_mRNA = result_mRNA[(result_mRNA$adj.P.Val < adj.P.Val & (result_mRNA$logFC>logFoldChange | result_mRNA$logFC<(-logFoldChange))),]
write.table(diffSig_mRNA,file="diffSig_logFC1_Padj0.01.txt",sep="\t",quote=F)

