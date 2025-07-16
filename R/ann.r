setwd
library(limma)
library(pheatmap)
library(ggplot2)
library(data.table)
library(dplyr)

data_GEO <- read.table("GSE66957_series_matrix.txt",header = T,check.names = F,sep = "\t")
GPL <- read.table("GPL15048.txt",header = T,check.names = F,sep = "\t")

colnames(GPL)[1] <- "ID_REF" 
data_GEO <- merge(GPL,data_GEO,by="ID_REF") 
data_GEO <- data_GEO[,-1]

rt=as.matrix(data_GEO) 
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)] 
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data) 
data_GEO <- data 

range(data)
write.csv(data,file="GSE66957.csv",row.names = T)


