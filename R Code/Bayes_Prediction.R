### Clear Console ###
cat("\014")

### Clear Environment ### 
rm(list = ls(all = TRUE))

### Load in the R libraries ###
library(BGLR)
library(MASS)
library(doParallel)
library(R.matlab)

### R Cancer/Genomics Packages From BioConductor ###
source("http://bioconductor.org/biocLite.R")

######################################################################################
######################################################################################
######################################################################################

### Run the MATLAB Code that creates the SECT Matrices for each image ###
system("matlab -nodisplay -r \"run('/home/lac55/GBM_Scripts/EC3D.m'); exit\"")

### Load in the structural array that holds the SECT Matrices ###
MRIs = readMat('/data/mukherjeelab/GBM/MRI_SECTs.mat')
SECTs = matrix(unlist(MRIs$MRIs[seq(2,length(MRIs$MRIs),2)]),nrow = length(MRIs$MRIs)/2,byrow = TRUE)
rownames(SECTs) = unlist(MRIs$MRIs[seq(1,length(MRIs$MRIs),2)])
SECTs = SECTs[,-c(101*(1:72),101*(0:72)+1)] #Get rid of noise

### Load in the GBM mRNA Gene Expression Data ###
library(diggitdata)
data(gbm.expression)
X_TCGA = exprs(gbmExprs)
colnames(X_TCGA) = substr(colnames(X_TCGA),start=1,stop=12)
G = t(X_TCGA[,which(colnames(X_TCGA)%in%rownames(SECTs))])
SECTs = SECTs[which(rownames(SECTs)%in%rownames(G)),]
G = G[match(rownames(G),rownames(SECTs)),]

### Load in the Morphometric and Volumetric Features ###
Morph = read.table("/data/mukherjeelab/GBM/Morphometric_Features.txt",header = TRUE)
rownames(Morph) = paste(Morph$Feature,Morph$Statistics,sep = "_")
colnames(Morph) = gsub("[.]", "-", as.character(colnames(Morph)))
Morph = t(Morph[,-c(1:2)])
Morph = Morph[which(rownames(Morph)%in%rownames(SECTs)),]

Geo = read.csv("/data/mukherjeelab/GBM/Volumetric_Features.csv")
rownames(Geo) = Geo$Patient.ID; Geo = Geo[,-1]
Geo = Geo[which(rownames(Geo)%in%rownames(SECTs)),]

SECTs = SECTs[which(rownames(SECTs)%in%rownames(Morph)),]
G = G[which(rownames(G)%in%rownames(Morph)),]
Geo = Geo[which(rownames(Geo)%in%rownames(Morph)),]

SECTs = SECTs[which(rownames(SECTs)%in%rownames(Geo)),]
G = G[which(rownames(G)%in%rownames(Geo)),]
Morph = Morph[which(rownames(Morph)%in%rownames(Geo)),]

G = G[match(rownames(G),rownames(SECTs)),]
Geo = Geo[match(rownames(Geo),rownames(SECTs)),]
Morph = Morph[match(rownames(Morph),rownames(SECTs)),]

dim(SECTs); dim(G); dim(Geo); dim(Morph);

######################################################################################
######################################################################################
######################################################################################

### Read in the Phenotypes ###
Phenos = read.csv("/data/mukherjeelab/GBM/TCGA_Clinical_Data.csv")

### Keep only disease free survival (DFS) and overall survival (OS) ###
Y = Phenos; rownames(Y) = as.character(Y$Patient.ID); Y = Y[,-1]
Y = Y[which(rownames(Y)%in%rownames(SECTs)),c(17,18)] 

### Check the Dimensionalities of the Data ###
dim(SECTs); dim(G); dim(Geo); dim(Morph); dim(Y)

######################################################################################
######################################################################################
######################################################################################

### Set up the Prediction (Classification) Parameters ###
setwd("/data/mukherjeelab/GBM/")

### Set a Random Seed for Reproducibility ###
set.seed(11151990)

#Define the percentage of data to use for splits
nsplit = 0.8

#Define the number of datasets to simulate (i.e. times to split)
ndatasets = 500
iter = 2e3
burn = 1e3

### Call the available cores accross the computer ###
registerDoParallel(cores = detectCores())

### Set up List ###
FMR = list()

for(j in 1:ncol(Y)){
  
  ### Select the Phenotype to model ###
  X = scale(G[!is.na(Y[,j]),],scale = FALSE)
  E = scale(SECTs[!is.na(Y[,j]),],scale = FALSE)
  V = scale(Geo[!is.na(Y[,j]),],scale = FALSE)
  M = scale(Morph[!is.na(Y[,j]),],scale = FALSE)
  y = scale(Y[!is.na(Y[,j]),j])
  
  KX = tcrossprod(X)/ncol(X)
  KE = tcrossprod(E)/ncol(E)
  KV = tcrossprod(V)/ncol(V)
  KM = tcrossprod(M)/ncol(M)
  
  ### Set up a list to save results ###
  Results = foreach(i = 1:ndatasets)%dopar%{
    
    ### Create the training and test sets ###
    ind = sample(1:length(y), size=nsplit*length(y), replace=FALSE)
    
    ### Subset the training and test sets ###
    yNA = y; yNA[-ind] = NA 
    
    ######################################################################################
    
    ### Model 1: Gene Expression Only ###
    ETA = list(list(K = KX, model="RKHS"))
    
    ### Run the Gibbs Sampler ###
    reg1 = BGLR(y=yNA, ETA=ETA, nIter=iter, burnIn=burn, verbose=FALSE)
    
    ### Get the posterior of the missing variables ###
    pred1 = reg1$yHat[-ind]
    
    ### Get Diagnostics ###
    RMSPE1 = sqrt(mean((y[-ind]-pred1)^2))
    COR1 = cor(y[-ind],pred1)
    
    ######################################################################################
    
    ### Model 2: Gene Expression + Morphometric Features ###
    ETA = list(list(K = KX, model="RKHS"),list(K = KM, model="RKHS"))
    
    ### Run the Gibbs Sampler ###
    reg2 = BGLR(y=yNA, ETA=ETA, nIter=iter, burnIn=burn, verbose=FALSE)
    
    ### Get the posterior of the missing variables ###
    pred2 = reg2$yHat[-ind]
    
    ### Get Diagnostics ###
    RMSPE2 = sqrt(mean((y[-ind]-pred2)^2))
    COR2 = cor(y[-ind],pred2)
    
    ######################################################################################
    
    ### Model 3: Gene Expression + Volumetric Features ###
    ETA = list(list(K = KX, model="RKHS"),list(K = KV, model="RKHS"))
    
    ### Run the Gibbs Sampler ###
    reg3 = BGLR(y=yNA, ETA=ETA, nIter=iter, burnIn=burn, verbose=FALSE)
    
    ### Get the posterior of the missing variables ###
    pred3 = reg3$yHat[-ind]
    
    ### Get Diagnostics ###
    RMSPE3 = sqrt(mean((y[-ind]-pred3)^2))
    COR3 = cor(y[-ind],pred3)
    
    ######################################################################################
    
    ### Model 4: Gene Expression + Volumetric Features + Morphometric Features ###
    ETA = list(list(K = KX, model="RKHS"),list(K = KM, model="RKHS"),list(K = KV, model="RKHS"))
    
    ### Run the Gibbs Sampler ###
    reg4 = BGLR(y=yNA, ETA=ETA, nIter=iter, burnIn=burn, verbose=FALSE)
    
    ### Get the posterior of the missing variables ###
    pred4 = reg4$yHat[-ind]
    
    ### Get Diagnostics ###
    RMSPE4 = sqrt(mean((y[-ind]-pred4)^2))
    COR4 = cor(y[-ind],pred4)
    
    ######################################################################################
    
    ### Model 5: Morphometric Features Only ###
    ETA = list(list(K = KM, model="RKHS"))
    
    ### Run the Gibbs Sampler ###
    reg5 = BGLR(y=yNA, ETA=ETA, nIter=iter, burnIn=burn, verbose=FALSE)
    
    ### Get the posterior of the missing variables ###
    pred5 = reg5$yHat[-ind]
    
    ### Get Diagnostics ###
    RMSPE5 = sqrt(mean((y[-ind]-pred5)^2))
    COR5 = cor(y[-ind],pred5)
    
    ######################################################################################
    
    ### Model 6: Volumetric Features Only ###
    ETA = list(list(K = KV, model="RKHS"))
    
    ### Run the Gibbs Sampler ###
    reg6 = BGLR(y=yNA, ETA=ETA, nIter=iter, burnIn=burn, verbose=FALSE)
    
    ### Get the posterior of the missing variables ###
    pred6 = reg6$yHat[-ind]
    
    ### Get Diagnostics ###
    RMSPE6 = sqrt(mean((y[-ind]-pred6)^2))
    COR6 = cor(y[-ind],pred6)
    
    ######################################################################################
    
    ### Model 7: Morphometric Features + Volumetric Features ###
    ETA = list(list(K = KV, model="RKHS"),list(K = KM, model="RKHS"))
    
    ### Run the Gibbs Sampler ###
    reg7 = BGLR(y=yNA, ETA=ETA, nIter=iter, burnIn=burn, verbose=FALSE)
    
    ### Get the posterior of the missing variables ###
    pred7 = reg7$yHat[-ind]
    
    ### Get Diagnostics ###
    RMSPE7 = sqrt(mean((y[-ind]-pred7)^2))
    COR7 = cor(y[-ind],pred7)
    
    ######################################################################################
    
    ### Model 8: SECT Summaries Only ###
    ETA = list(list(K = KE, model="RKHS"))
    
    ### Run the Gibbs Sampler ###
    reg8 = BGLR(y=yNA, ETA=ETA, nIter=iter, burnIn=burn, verbose=FALSE)
    
    ### Get the posterior of the missing variables ###
    pred8 = reg8$yHat[-ind]
    
    ### Get Diagnostics ###
    RMSPE8 = sqrt(mean((y[-ind]-pred8)^2))
    COR8 = cor(y[-ind],pred8)
    
    ######################################################################################
    
    ### Model 9: Gene Expression + SECT Summaries ###
    ETA = list(list(K = KX, model="RKHS"),list(K = KE, model="RKHS"))
    
    ### Run the Gibbs Sampler ###
    reg9 = BGLR(y=yNA, ETA=ETA, nIter=iter, burnIn=burn, verbose=FALSE)
    
    ### Get the posterior of the missing variables ###
    pred9 = reg9$yHat[-ind]
    
    ### Get Diagnostics ###
    RMSPE9 = sqrt(mean((y[-ind]-pred9)^2))
    COR9 = cor(y[-ind],pred9)
    
    ######################################################################################
    
    ### Model 10: SECT Summaries + Morphometric Features ###
    ETA = list(list(K = KE, model="RKHS"),list(K = KM, model="RKHS"))
    
    ### Run the Gibbs Sampler ###
    reg10 = BGLR(y=yNA, ETA=ETA, nIter=iter, burnIn=burn, verbose=FALSE)
    
    ### Get the posterior of the missing variables ###
    pred10 = reg10$yHat[-ind]
    
    ### Get Diagnostics ###
    RMSPE10 = sqrt(mean((y[-ind]-pred10)^2))
    COR10 = cor(y[-ind],pred10)
    
    ######################################################################################
    
    ### Model 11: SECT Summaries + Volumetric Features ###
    ETA = list(list(K = KE, model="RKHS"),list(K = KV, model="RKHS"))
    
    ### Run the Gibbs Sampler ###
    reg11 = BGLR(y=yNA, ETA=ETA, nIter=iter, burnIn=burn, verbose=FALSE)
    
    ### Get the posterior of the missing variables ###
    pred11 = reg11$yHat[-ind]
    
    ### Get Diagnostics ###
    RMSPE11 = sqrt(mean((y[-ind]-pred11)^2))
    COR11 = cor(y[-ind],pred11)
    
    ######################################################################################
    
    ### Model 12: SECT Summaries + Morphometric Features + Volumetric Features ###
    ETA = list(list(K = KE, model="RKHS"),list(K = KM, model="RKHS"),list(K = KV, model="RKHS"))
    
    ### Run the Gibbs Sampler ###
    reg12 = BGLR(y=yNA, ETA=ETA, nIter=iter, burnIn=burn, verbose=FALSE)
    
    ### Get the posterior of the missing variables ###
    pred12 = reg12$yHat[-ind]
    
    ### Get Diagnostics ###
    RMSPE12 = sqrt(mean((y[-ind]-pred12)^2))
    COR12 = cor(y[-ind],pred12)
    
    ######################################################################################
    
    ### Save The Results ###
    c(RMSPE1,RMSPE2,RMSPE3,RMSPE4,RMSPE5,RMSPE6,RMSPE7,RMSPE8,RMSPE9,RMSPE10,RMSPE11,RMSPE12,COR1,COR2,COR3,COR4,COR5,COR6,COR7,COR8,COR9,COR10,COR11,COR12)
  }
  
  ### Create a matrix to save results ###
  Final = matrix(unlist(Results),nrow = ndatasets,ncol = 24,byrow = TRUE)
  colnames(Final) = c(paste("RMSE",1:12,sep="_"),paste("COR",1:12,sep="_"))
  
  FMR[[j]] = Final
  
  cat("Completed Phenotype",j,"\n")
}

#colMeans(FMR[[1]])
#boxplot(FMR[[1]][,2:5]-FMR[[1]][,1])

### Save Results ###
save(FMR,file = "Pred_Results_Bayes.RData")