#' @name computeConditionalSimilarities
#' @title Compute conditional PPI similarities
#'
#' @description Score  similarity of interacting protein CF/MS profiles in each
#' condition.
#' 
#'  There are two scores computed:
#' * ***Co-elution***: Pearson correlation of CF/MS profiles
#'
#' * ***Co-abundance***:  Log2 fold-change of summed CF/MS profile
#'  intensities. Since CF/MS-derived PPI are undirected, fold-change is
#'  fixed as ratio between lower intensity and higher one of the two 
#'  interacting proteins.
#'    
#' @param proteins List of lists. One list per replicate.
#'  Each replicate list consists of one protein-by-fraction data frame
#'  per condition.
#'  
#' @param conditions Character vector with names of experimental conditions.
#' 
#' @param referencePPI Character data frame with two columns corresponding to
#' the protein IDs of the interacting proteins in the reference interactome. 
#' Protein IDs must match the format of those in the \code{proteins} input.
#' 
#' @param logTransformed. Logical. Whether the input protein profiles were
#' log2-transformed. Default TRUE.
#'
#' @return Data frame. Similarity score for each PPI in each condition.
#' 
#' @examples
#' conditionalSimilarities <- computeConditionalSimilarities(proteins,
#'  conditions,
#'  referencePPI)
#'
#' @import dplyr magrittr WGCNA reshape2
NULL
#' @export
computeConditionalSimilarities <- function(proteins,
                                           conditions,
                                           referencePPI,
                                           logTransformed = TRUE){
  
  #compute co-elution scores
  coelution <- lapply(conditions, function(curCond){
    #iterate over replicates
    conditionalSimilarities <- lapply(proteins, function(curRep){
      #get profiles of this condition-replicate
      curProfiles <- curRep[[curCond]]
      
      #compute Pearson correlation between all pairs of proteins in this condition
      allCors <- WGCNA::cor(t(curProfiles))
      
      #get distances of reference protein pairs
      refCors <- diag(allCors[referencePPI[,1], referencePPI[,2]])
      
      return(refCors)
    }) %>% 
      #average over replicates
      do.call("cbind",.) %>% 
      rowMeans()
    
    return(conditionalSimilarities)
  }) %>% 
    #add names of conditions
    setNames(conditions) %>% 
    as.data.frame() %>% 
    #add interacting protein IDs
    cbind(referencePPI, .)
  
  #compute co-abundance scores
  coabundance <- lapply(conditions, function(curCond){
    #iterate over replicates
    conditionalFoldChanges <- lapply(proteins, function(curRep){
      #get profiles of this condition-replicate
      curProfiles <- curRep[[curCond]]
      
      #reverse log2-transformation
      if(logTransformed)
        curProfiles <- (2^curProfiles) - 1
      
      #sum intensity across fractions for each protein
      intensitySums <- rowSums(curProfiles)
      
      #compute fold-changes between all pairs of proteins
      allFoldChanges <- outer(intensitySums, intensitySums, FUN = "/")
      
      #get distances of reference protein pairs
      refFoldChanges <- diag(allFoldChanges[referencePPI$V1, referencePPI$V2])
      
      #convert fold-changes so that it is always min/max
      refFoldChanges[refFoldChanges > 1] <- 1 / refFoldChanges[refFoldChanges > 1]
      
      #log2-transform
      refFoldChanges <- log2(refFoldChanges)
      
      return(refFoldChanges)
    }) %>% 
      #average over replicates
      do.call("cbind",.) %>% 
      rowMeans()
    
    return(conditionalFoldChanges)
  }) %>% 
    #add names of conditions
    setNames(conditions) %>% 
    as.data.frame() %>% 
    #add interacting protein IDs
    cbind(referencePPI, .)
  
  #combine similarity scores into one using PCA
  combinedSimilarities <- cbind(melt(ppiCors)$value, melt(foldChanges)$value) %>% 
    #scale columns
    scale() %>% 
    #PCA
    prcomp(.) %>% 
    #retain 1st PC only
    .$x %>% .[,1] %>% 
    #add PPI/condition info
    cbind(melt(ppiCors)[,1:3], .) %>% 
    set_colnames(c("V1", "V2", "condition", "score")) %>% 
    dcast(V1 + V2  ~ condition)
  
  return(combinedSimilarities)
}

