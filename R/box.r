setwd
library(reshape2)
library(ggpubr)
inputFile="input.txt"      
outFile="boxplot.pdf"      

rt=read.table(inputFile,sep="\t",header=T,check.names=F,row.names=1)
x=colnames(rt)[1]
colnames(rt)[1]="Type"

data=melt(rt,id.vars=c("Type"))
colnames(data)=c("Type","Gene","Expression")

p=ggboxplot(data, x="Gene", y="Expression", color = "Type", 
            ylab="",
            xlab="",
            legend.title=x,
            palette = c("blue","red"),
            width=0.6, add = "none")
p=p+rotate_x_text(60)
p1=p+stat_compare_means(aes(group=Type),
                        method="wilcox.test",
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif")

pdf(file=outFile, width=8, height=8)
print(p1)
dev.off()