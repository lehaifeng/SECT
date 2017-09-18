colMeans(FMR[[1]][,1:4]); apply(FMR[[1]][,1:4],2,var) 

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
round(cbind(colMeans(OS),apply(OS^2,2,sd)/sqrt(n.splits)),3)

######################################################################################
######################################################################################
######################################################################################

round(colMeans(FMR[[2]][,5:8]),3)
round(colMeans(FMR[[1]][,5:8]),3)
