### Clear Console ###
cat("\014")

### Clear Environment ### 
rm(list = ls(all = TRUE))

### Load in the R libraries ###
library(coda)
library(MASS)
library(doParallel)
library(Rcpp)
library(RcppArmadillo)
library(RcppParallel)
library(BGLR)
library(monomvn)
library(kernlab)
library(adegenet)
library(RColorBrewer)

######################################################################################
######################################################################################
######################################################################################

### Set the Working Directory ###
setwd("~/Dropbox/Desktop/SECT Revision/Analysis/")

### Load in the Prediction Results ###
#load("SECT_Gauss_Results.RData");
#load("SECT_Linear_Results.RData");
#load("SECT_Sigmoid_Results.RData");
load("SECT_Log_Results.RData");

### Set the Number of Splits Used in the Simulation ###
n.splits = nrow(FMR[[1]])

DFS = FMR[[2]][,1:4]
OS = FMR[[1]][,1:4]

countMSE = matrix(0,nrow = n.splits,ncol = 4)
for(i in 1:nrow(FMR[[1]])){
  countMSE[i,which(DFS[i,] == min(DFS[i,]))] = 1
}

colMeans(countMSE)

countMSE = matrix(0,nrow = n.splits,ncol = 4)
for(i in 1:nrow(FMR[[1]])){
  countMSE[i,which(OS[i,] == min(OS[i,]))] = 1
}

colMeans(countMSE)

round(cbind(colMeans(DFS),apply(DFS,2,sd)/sqrt(n.splits)),3)
round(cbind(colMeans(OS),apply(OS,2,sd)/sqrt(n.splits)),3)

######################################################################################
######################################################################################
######################################################################################

### Create Method Groups ###
models = c("Expression","Morphometric","Volume","SECT")
groups = c()
for(i in 1:length(models)){
  groups = factor(c(groups,rep(models[i],100)))
}

x = factor(1:length(models))
col.info = any2col(x,col.pal = funky)

### Create Boxplots of Relative
pdf("Related_Pred_Fig.pdf",width=11, height=4.25)
par(mar=c(5.1, 4.1, 5, 2), xpd=TRUE)
bp=boxplot(as.numeric(DFS)~groups,col = col.info$col,frame = FALSE, ylim = c(0,2.5),xlim = c(0,10),names = FALSE, at = 1:4,xaxt = "n",ylab = "MSPE",xlab = "Scenarios",boxwex=0.8,medlwd = 1.5,outpch = 1)
boxplot(as.numeric(OS)~groups,col = col.info$col,frame = FALSE,names = FALSE, at = 6:9, add = TRUE,xaxt = "n",boxwex=0.8,medlwd = 1.5,outpch = 1)
axis(side = 1,at=c(0,5,15,25,30.5),labels = c("","I","II","III",""))
segments(x0 = 10, x1 = 10,y0 = 0.3, y1 =4, lty = 2, col = "grey", lwd = 1.5)
L = legend(x = 'top', legend = c("BRR","BL","BLMM","SVM",expression("BAKR ("*h==5*")"),expression("BAKR ("*h==2*")"),expression("BAKR ("*h==1*")"),expression("BAKR ("*h==0.5*")"),expression("BAKR ("*h==0.01*")")), col=col.info$col, ncol=5, bty='n', pch=15, inset=c(0,-0.3),cex = 0.8)
dev.off()
