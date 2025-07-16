setwd
library(glmnet)
library(tidyr)
lasso.v <- read.table("input.txt",header=T,sep="\t",row.names=1)
group.v <- c(rep("case",57),rep('control',12)) %>% factor(.,levels = c("case","control"),ordered = F)
design.v <- model.matrix(~factor(group.v))   
lasso.v<-cbind(design.v[,2],lasso.v)
colnames(lasso.v)[1]<- 'group'
lasso.v <- as.matrix(lasso.v)
set.seed(123)
fit_stage.v <-cv.glmnet(lasso.v[,-1],lasso.v[,'group'],
                        family="binomial",type.measure="default",nfolds = 3, gamma = 1)
fit.v <- fit_stage.v$glmnet.fit
pdf(file = "lasso.Binomial.Deviance.pdf",height = 8,width = 8)
par(mgp = c(2.5,1,0),mar=c(4,5,3,3))
plot(fit_stage.v,xlab='Log Lambda',cex.lab = 1)+
  text(x = log(fit_stage.v$lambda.min),y = 0.8,
       paste('Lambda.min\n',round(fit_stage.v$lambda.min,4)),cex=0.8,adj=0.1)+
  text(x = log(fit_stage.v$lambda.1se),y = 1,
       paste('Lambda.lse\n',round(fit_stage.v$lambda.1se,4)),cex=0.8)
dev.off()

pdf(file = "lasso.voefficients.venalty.pdf",height = 8,width = 8)
par(mgp = c(2.5,1,0),mai=c(1,1,0.5,0.5))
plot(fit.v, xvar="lambda",cex.lab = 1)+
  abline(v = c(log(fit_stage.v$lambda.min), log(fit_stage.v$lambda.1se)),lty=2)+
  text(x = log(fit_stage.v$lambda.min),y = 0,
       paste('Lambda.min\n',round(fit_stage.v$lambda.min,4)),cex=0.8,adj=0.5)+
  text(x = log(fit_stage.v$lambda.1se),y = 2,
       paste('Lambda.lse\n',round(fit_stage.v$lambda.1se,4)),cex=0.8,adj=0.5)
dev.off()

lam = fit_stage.v$lambda.min
pre_train<-predict(fit_stage.v,type = "class", newx=lasso.v[,-1], s=lam)
acc_table_tr<-table(pre_train,lasso.v[,'group'])
colnames(acc_table_tr)<-c("normal","case")
rownames(acc_table_tr)<-c("normal","case")
suppressMessages(library(kableExtra))
kable(acc_table_tr, "html",caption = '<center>**表. LASSO(Training set) **</center>') %>%
  kable_styling(bootstrap_options = "striped", full_width = F)

coefficients<-coef(fit_stage.v,s=lam)
Active.Index<-coefficients@i
Active.voefficients<-coefficients@x
Active.name.train <- colnames(lasso.v[,-1])[Active.Index[-1]]
write.csv(Active.name.train,"lasso.csv")