## Script for computing pairwise protein distances across conditions ##

#required packages
library(pbapply) #parallelize operations
library(magrittr) #data wrangling
library(dplyr) #data wrangling
library(purrr) #data wrangling

#input: List of conditional measurements. One list per replicate. In each replicate's list: one protein-by-fraction table per condition.
#input: Distance metric for computing distances between proteins.
#input: Whether there are replicate measurements to account for. Default TRUE.
#input: Number of cores for parallelization of operations. Default 6.
#output: Matrix with distance between protein pairs in each condition. One column per condition.
computePPI <- function(cofracTables, metric = c('wasserstein', 'euclidean', 'pearson'), replicates = TRUE, cores = 6){
  #names of experimental conditions
  conditions <- names(cofracTables)
  
  #merge replicate measurements into one table for each condition
  if(replicates == TRUE){
    #names of experimental conditions
    conditions <- names(cofracTables[[1]])
    
    cofracTables <- lapply(conditions, function(curCond){
      curReps <- lapply(cofracTables, `[[`, curCond)
      #remove zero-intensity profiles
      curReps <- lapply(curReps, function(x) return(x[rowSums(x[,-1]) > 0,]))
      #retain reproducible proteins
      reproducibleProts <-  Reduce(intersect, lapply(curReps, "[", ,1))
      curReps <- do.call("rbind", curReps)
      curReps <- curReps[curReps[,1] %in% reproducibleProts,]
      
      return(curReps)
    })
    names(cofracTables) <- conditions
  }
  
  #retain only the proteins common to all experimental conditions
  toRetain <- Reduce(intersect, lapply(cofracTables, "[", ,1))
  cofracTables <- lapply(cofracTables, function(curCond){
    return(curCond[curCond[,1] %in% toRetain,])
  })
  
  ## Compute pairwise distances ##
  print('Computing distances between all pairs of proteins in each condition...')
  ppi <- pblapply(conditions, function(curCond){
    #get relevant measurements
    curProts <- cofracTables[[curCond]]
    
    #compute distances
    ppi <- computeDistances(curProts, metric)
    
    #average distances for each unique pairing (relevant if there are replicates)
    if(replicates == TRUE){
      ppi <- ppi %>%
        group_by(id1, id2) %>%
        dplyr::summarise(distance = mean(distance))
      
      #divide by noise
      noise <- ppi[ppi$id1==ppi$id2,]
      noise <- setNames(unlist(noise$distance), noise$id1)
      ppi$noise <- (noise[ppi$id1] + noise[ppi$id2]) / 2
      ppi$distance <- ppi$distance / ppi$noise
      
      colnames(ppi)[1:3] <- c("id1", "id2", curCond)
      
      return(ppi[,-4])
    }
    else{
      colnames(ppi) <- c("id1", "id2", curCond)
      return(ppi)
    }
  }, cl = cores)
  names(ppi) <- conditions
  
  #merge conditional scores into one table
  ppi <-  ppi %>% 
    purrr::reduce(left_join, by = c("id1","id2"))
  
  return(ppi)
}
