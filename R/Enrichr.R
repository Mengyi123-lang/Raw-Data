setwd
library('enrichR')
library(clusterProfiler)
library(org.Hs.eg.db)
#dbs <- listEnrichrDbs()
dbs<- c("DSigDB")
symbol=read.csv("input.csv",header = F)
enrichr<- enrichr(symbol$V1, dbs)
write.csv(enrichr$DSigDB,"EnrichR.csv")