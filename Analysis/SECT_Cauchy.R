### Clear Console ###
cat("\014")

### Clear Environment ### 
rm(list = ls(all = TRUE))

### Load in the R libraries ###
library(BGLR)
library(doParallel)
library(gdata)
library(diggitdata)
library(Rcpp)
library(RcppArmadillo)
library(RcppParallel)

### Load in the C++ BAKR functions ###
sourceCpp("/home/lcrawfo1/Scripts/BAKR/BAKRGibbs.cpp")

### R Cancer/Genomics Packages From BioConductor ###
source("http://bioconductor.org/biocLite.R")

######################################################################################
######################################################################################
######################################################################################

### Run the MATLAB Code that creates the EC Matrices for each image ###
#system("/Users/lorincrawford/Applications/MATLAB_R2017a.app/bin/matlab -nodisplay -r \"run('~/Dropbox/Columbia Radiogenomics/Software/CompEC.m'); exit\"")

### Load in the structural array that holds the EC Matrices ###
#Shapes = readMat('~/Dropbox/Columbia Radiogenomics/Data/ECs.mat')
#ECs = matrix(unlist(Shapes$Shapes[seq(2,length(Shapes$Shapes),2)]),nrow = length(Shapes$Shapes)/2,byrow = TRUE)

### Run the R Code that creates the EC Matrices for each image ###
#source("~/Dropbox/Columbia Radiogenomics/Software/EC3D.R")
#setwd("~/Dropbox/Columbia Radiogenomics/Data")
#startdir = getwd()
#in.dir = "MITKSegmentations"
#out.file = "~/Dropbox/Columbia Radiogenomics/Data/MRIECs.RData"
#ecf = ecf(in.dir = dir,out.file = out.file,img.dir = "baseline/Segmentations/enh",first.only = FALSE)

### Load in the structural array that holds the EC Matrices ###
load('/home/lcrawfo1/Data/GBM/MRIECs.RData')
nrot = ncol(MRI_list[[1]]$EC); stepsize = nrow(MRI_list[[1]]$EC)
ECs = matrix(nrow = length(MRI_list),ncol = nrot*stepsize)
rownames(ECs) = 1:nrow(ECs)
dim(ECs)

for(i in 1:nrow(ECs)){
  ECs[i,] = c(MRI_list[[i]]$EC)
  rownames(ECs)[i] = MRI_list[[i]]$name
}
ECs = ECs[,-c(101*(1:nrot),101*(0:nrot)+1)]  
ECs[1:10,1:10]

library(diggitdata)
data(gbm.expression)
X_TCGA = exprs(gbmExprs)
colnames(X_TCGA) = substr(colnames(X_TCGA),start=1,stop=12)
G = t(X_TCGA[,which(colnames(X_TCGA)%in%rownames(ECs))])
ECs = ECs[which(rownames(ECs)%in%rownames(G)),]
G = G[match(rownames(ECs),rownames(G)),]

### Load in the Morphometric and Volumetric Features ###
Morph = read.table("/home/lcrawfo1/Data/GBM/patient_statistics.txt",header = TRUE)
rownames(Morph) = paste(Morph$Feature,Morph$Statistics,sep = "_")
colnames(Morph) = gsub("[.]", "-", as.character(colnames(Morph)))
Morph = t(Morph[,-c(1:2)])
Morph = Morph[which(rownames(Morph)%in%rownames(ECs)),]

Geo = read.xls("/home/lcrawfo1/Data/GBM/TCGA-geometric-REDO.xls")
rownames(Geo) = Geo$Patient.ID; Geo = as.matrix(Geo[,-c(1,6:11)])
Geo = Geo[which(rownames(Geo)%in%rownames(ECs)),]

ECs = ECs[which(rownames(ECs)%in%rownames(Morph)),]
G = G[which(rownames(G)%in%rownames(Morph)),]
Geo = Geo[which(rownames(Geo)%in%rownames(Morph)),]

ECs = ECs[which(rownames(ECs)%in%rownames(Geo)),]
G = G[which(rownames(G)%in%rownames(Geo)),]
Morph = Morph[which(rownames(Morph)%in%rownames(Geo)),]

G = G[match(rownames(ECs),rownames(G)),]
Geo = Geo[match(rownames(ECs),rownames(Geo)),]
Morph = Morph[match(rownames(ECs),rownames(Morph)),]

dim(ECs); dim(G); dim(Geo); dim(Morph);

######################################################################################
######################################################################################
######################################################################################

### Read in the Phenotypes ###
Phenos = read.csv("/home/lcrawfo1/Data/GBM/TCGA genomic data.csv")

### Denote the Two classes as 1s or 0s Depending on the class ###
Y = Phenos; rownames(Y) = as.character(Y$Patient.ID); Y = Y[,-1]
Y = Y[which(rownames(Y)%in%rownames(ECs)),17:18]

ECs = ECs[which(rownames(ECs)%in%rownames(Y)),]
G = G[which(rownames(G)%in%rownames(Y)),]
Morph = Morph[which(rownames(Morph)%in%rownames(Y)),]
Geo = Geo[which(rownames(Geo)%in%rownames(Y)),]

### NOTE: [[Odd]] has the shape names and [[Even]] has the EC Matrices ###

### Check the Dimensionalities of the Data ###
dim(ECs); dim(G); dim(Geo); dim(Morph); dim(Y)

######################################################################################
######################################################################################
######################################################################################

### Models to Consider: (1) Bayesian Lasso; (2) Bayesian Ridge; (3) Bayes B; (4) Bayes Cpi; (5) Bayesian BLUP###

### Set the working directory for results ###
setwd("/home/lcrawfo1/Results/SECT")

### Set the random seed for reproducibility ###
set.seed(11151990)

#Define the percentage of data to use for splits
nsplit = 0.75

#Define the number of datasets to simulate (i.e. times to split)
ndatasets = 100
iter = 2e4
burn = 1e4
thin = 10
h=1

### Call the available cores accross the computer ###
registerDoParallel(cores = detectCores())

### Set up List ###
FMR = list()

######################################################################################
######################################################################################
######################################################################################

### Run the Analysis ###
for(j in 1:ncol(Y)){
  
  ### Choose Phenotype ###
  y = scale(Y[!is.na(Y[,j]),j])
  
  ### Set up a list to save results ###
  Results = foreach(i = 1:ndatasets)%dopar%{
    
    ### Create the training and test sets ###
    ind = sample(1:length(y), size=nsplit*length(y), replace=FALSE)
    
    ### Subset the training and test sets ###
    yNA = y; yNA[-ind] = NA
    
    ######################################
    
    X = Geo[!is.na(Y[,j]),]
    
    K = tanh(LinearKernel(t(X),h=ncol(X)))
    
    ### Center and Scale K_tilde ###
    n=nrow(K)
    v=matrix(1, n, 1)
    M=diag(n)-v%*%t(v)/n
    Kn=M%*%K%*%M
    Kn=Kn/mean(diag(Kn))
    
    ### Bayesian Kernel Ridge Regression ###
    ETA = list(list(K=Kn, model="RKHS"))
    mod = BGLR(y=yNA, ETA=ETA, nIter=iter, burnIn=burn, verbose=FALSE,thin = thin)
    
    RMSPE.Geo = sqrt(mean((y[-ind]-mod$yHat[-ind])^2))
    R2_Geo = 1-(sum((y[ind]-mod$yHat[ind])^2)/sum((y[ind]-mean(y[ind]))^2))
    
    ######################################################################################
    ######################################################################################
    ######################################################################################
    
    ### Choose Phenotype ###
    X = Morph[!is.na(Y[,j]),]
    
    K = tanh(LinearKernel(t(X),h=ncol(X)))
    
    ### Center and Scale K_tilde ###
    n=nrow(K)
    v=matrix(1, n, 1)
    M=diag(n)-v%*%t(v)/n
    Kn=M%*%K%*%M
    Kn=Kn/mean(diag(Kn))
    
    ### Bayesian Kernel Ridge Regression ###
    ETA = list(list(K=Kn, model="RKHS"))
    mod = BGLR(y=yNA, ETA=ETA, nIter=iter, burnIn=burn, verbose=FALSE,thin = thin)
    
    RMSPE.Morph = sqrt(mean((y[-ind]-mod$yHat[-ind])^2))
    R2_Morph = 1-(sum((y[ind]-mod$yHat[ind])^2)/sum((y[ind]-mean(y[ind]))^2))
    
    ######################################################################################
    ######################################################################################
    ######################################################################################
    
    X = G[!is.na(Y[,j]),]
    
    K = tanh(LinearKernel(t(X),h=ncol(X)))
    
    ### Center and Scale K_tilde ###
    n=nrow(K)
    v=matrix(1, n, 1)
    M=diag(n)-v%*%t(v)/n
    Kn=M%*%K%*%M
    Kn=Kn/mean(diag(Kn))
    
    ### Bayesian Kernel Ridge Regression ###
    ETA = list(list(K=Kn, model="RKHS"))
    mod = BGLR(y=yNA, ETA=ETA, nIter=iter, burnIn=burn, verbose=FALSE,thin = thin)
    
    RMSPE.G = sqrt(mean((y[-ind]-mod$yHat[-ind])^2))
    R2_G = 1-(sum((y[ind]-mod$yHat[ind])^2)/sum((y[ind]-mean(y[ind]))^2))
    
    ######################################################################################
    ######################################################################################
    ######################################################################################
    
    ### Choose Phenotype ###
    X = ECs[!is.na(Y[,j]),]; X = scale(X)
    
    K = tanh(LinearKernel(t(X),h=ncol(X)))
    
    ### Center and Scale K_tilde ###
    n=nrow(K)
    v=matrix(1, n, 1)
    M=diag(n)-v%*%t(v)/n
    Kn=M%*%K%*%M
    Kn=Kn/mean(diag(Kn))
    
    ### Bayesian Kernel Ridge Regression ###
    ETA = list(list(K=Kn, model="RKHS"))
    mod = BGLR(y=yNA, ETA=ETA, nIter=iter, burnIn=burn, verbose=FALSE,thin = thin)
    
    RMSPE.ECs = sqrt(mean((y[-ind]-mod$yHat[-ind])^2))
    R2_ECs = 1-(sum((y[ind]-mod$yHat[ind])^2)/sum((y[ind]-mean(y[ind]))^2))
    
    ######################################
    
    ### Save The Results ###
    c(RMSPE.G,RMSPE.Geo,RMSPE.Morph,RMSPE.ECs,R2_G,R2_Geo,R2_Morph,R2_ECs)
  }
  ### Create a matrix to save results ###
  Final = matrix(unlist(Results),nrow = ndatasets,ncol = 8,byrow = TRUE)
  mod.names = c("Expression","Geo","Morph","ECs")
  colnames(Final) = c(paste("RMSEP",mod.names,sep="_"),paste("R2",mod.names,sep="_"))
  
  FMR[[j]] = Final
  
  cat("Completed Phenotype",j,"\n")
}

colMeans(FMR[[1]][,1:4]); 
colMeans(FMR[[2]][,1:4])

### Save Results ###
save(FMR,file = "SECT_Sigmoid_Results.RData")

