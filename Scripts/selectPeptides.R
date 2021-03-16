## Script for selecting representative peptides for each protein ##

#required packages
library(pbapply) #parallelize operations
library(magrittr) #data wrangling
library(dplyr) #data wrangling
library(data.table) #data wrangling
library(igraph) #clustering peptides

## Helper function: compute distances between sibling peptides ##
# input: profile-by-fraction matrix, first column has peptide IDs, second column has protein IDs and must be named 'protein'
# output: list of sibling peptide distance matrices, one per protein
siblingDistances <- function(pepTable, metric = c('wasserstein', 'euclidean', 'pearson'), cores = 6){
  
  #split input table into tables of sibling peptides
  profiles <- split(pepTable, with(pepTable, protein), drop=TRUE)
  
  #compute all pairwise distances for each group of sibling peptides
  print("Computing distances between sibling peptides for filtering...")
  sibDistances <- pblapply(profiles, function(siblings){
    sibs <- computeDistances(cofracTable = peps$replicate_1$LB, metric = metric)
    return(sibs)
  }, cl = cores) 
  
  return(sibDistances)
}


## Filter peptides based on similarity to sibling peptides ##

#input: Profile-by-fraction matrix. First column has peptide IDs, second column has protein IDs and must be named 'protein'
#input: Distance metric for sibling peptides.
#input: Method for filtering peptides. One of 'quantile' or 'clustering'.
#input: Distance cutoff threshold as quantile (e.g. 0.95 would filter out the top 5% peptides with highest distance)
filterOutlierPeptides <- function(pepTable,
                                  metric = c('wasserstein', 'euclidean', 'pearson'),
                                  method = c("quantile", "clustering"),
                                  threshold = 0.95,
                                  cores = 6){
  #compute distances between sibling peptides
  sibDistances <- siblingDistances(pepTable, metric, cores)
  
  #filter outlier peptides
  if(method == "quantile"){
    
    #collapse sibling peptide distances into one table
    collapsedDis <- do.call("rbind", sibDistances)
    
    #in case of 'infinity' distance, set to max distance
    collapsedDis$distance[is.infinite(collapsedDis$distance)] <- max(collapsedDis$distance[is.finite(collapsedDis$distance)])
    
    #list of all unique peptides
    allPeps <- pepTable[,1]
    
    #get average distance for each peptide
    print("Computing average agreement of each peptide with its siblings...")
    sibAgreement <- pblapply(allPeps, function(pep){
      return(mean(unlist(collapsedDis[collapsedDis$id1==pep | collapsedDis$id2==pep, "distance"]), na.rm = TRUE))
    },cl = cores)
    names(sibAgreement) <- allPeps
    
    #retain peptides which fall below threshold
    sibAgreement[unlist(sibAgreement) > quantile(unlist(sibAgreement), threshold, na.rm = TRUE)] <- NULL
    sibAgreement[is.na(unlist(sibAgreement))] <- NULL
    
    #return filtered peptide table
    return(pepTable[pepTable[,1] %in% names(sibAgreement),])
  }
  
  if(method == "clustering"){
    print("Filtering outlier peptides using clustering...")
    #iterate over proteins and decide which peptides to retain
    clusteredPeps <- pblapply(sibDistances, function(siblings){
      #remove any NA values (non-reproducible peptides)
      siblings <- siblings[complete.cases(siblings),]
      
      #if single-peptide protein return null
      if(length(union(siblings$id1, siblings$id2)) <= 1){
        return(NULL)
      }
      
      #convert to distance matrix via igraph
      distMatrix <- igraph::graph.data.frame(siblings[,1:2], directed = FALSE)
      E(distMatrix)$weight <- unlist(siblings[,3])
      distMatrix <- as.dist(as.matrix(as_adjacency_matrix(distMatrix, attr = "weight", sparse = TRUE)))
      
      #average-linkage hierarchical clustering with 2 clusters
      pepsClust <- cutree(hclust(distMatrix, method = "average"), k = 2)
      
      #get peptides belonging to largest cluster
      toKeep <- names(pepsClust)[pepsClust == tail(names(sort(table(pepsClust))), 1)]
      
      #return peptides to keep
      return(toKeep)
      
    }, cl = cores)
    
    #keep peptides after filtering outliers
    return(pepTable[pepTable[,1] %in% unlist(clusteredPeps),])
  }
}


## Select top N peptides ##
#input: Profile-by-fraction matrix. First column has peptide IDs, second column has protein IDs and must be named 'protein'
#input: Number of top-intensity peptides to retain as representative of each protein 'n'. Default number is 3.
selectTopPeptides <- function(pepTable, n = 3){
  #sum up intensity of peptides
  pepSums <- rowSums(pepTable[,-c(1,2)])
  
  #create table with peptide, parent protein, and sum of intensity
  pepSums <- data.table(cbind(pepTable[,c(1,2)], pepSums), key = "protein")
  
  #in case of replicates, keep max sum across replicates for each peptide
  pepSums <- pepSums %>%
    group_by_at(1) %>% 
    filter(pepSums == max(pepSums))
  
  #find top N peptides for each protein  
  topPeps <- pepSums %>%
    arrange(desc(pepSums)) %>%
    group_by(protein) %>% 
    slice(1:n)
  
  #retain only the top N peptides for each protein
  filteredPeps <- pepTable[unlist(pepTable[,1]) %in% unlist(topPeps[,1]),]
  
  return(filteredPeps)
}
