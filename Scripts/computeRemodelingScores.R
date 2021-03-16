## Script for computing PPI remodeling scores across conditions ##

#required packages
library(pbapply) #parallelize operation

#input: Matrix with distance between protein pairs in each condition. First two columns are protein IDs, followed by one column per condition.
#input: Number of cores for parallelization of operations. Default 6.
#output: Matrix with remodeling scores for each pair of proteins and each pair of conditions.
computeRemodelingScores <- function(ppi, cores = 6){
  #list of condition names
  conditions <- colnames(ppi)[-c(1,2)]
  
  #generate matrix of all pairs of conditions
  condPairs <- t(combn(conditions, 2))
  
  #compute remodeling score as the difference in distances between each pair of conditions
  print("Computing remodeling scores between each pair of conditions...")
  remodelScores <- pbapply(condPairs, 1, function(pair){
    return(abs(ppi[,pair[1]] - ppi[,pair[2]]))
  }, cl = cores)
  
  #set column names
  colnames(remodelScores) <- apply(condPairs, 1, paste, collapse="_")
  
  #attach protein IDs
  remodelScores <- cbind(ppi[,1:2], remodelScores)
  
  return(remodelScores)
}
