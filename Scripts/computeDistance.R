## Compute distances between co-frac profiles using different metrics ##

#required packages
library(dplyr) #data wrangling
library(reshape2) #data wrangling
library(data.table) #data wrangling

#input: profile-by-fraction matrix, first column has IDs
#output: data frame with 3 columns: id 1, id 2, distance
computeDistances <- function(cofracTable, metric = c('wasserstein', 'euclidean', 'pearson')){
  #sort profile rows by ID
  cofracTable <- cofracTable[order(as.character(unlist(cofracTable[,1]))),]
  
  #store profile (protein/peptide) IDs in variable
  ids <- as.character(unlist(cofracTable[,1]))

  #keep intensity values only
  cofracTable <- dplyr::select_if(cofracTable, is.numeric)
  
  #wasserstein distance
  if(metric == "wasserstein"){
    # 1D Wasserstein distance is computed as the manhattan distance between the CDFs #
    
    #normalize profiles to sum to 1
    normTable <- t(apply(cofracTable, 1, function(profile){
      return(profile/sum(profile))
    }))
    
    #take running sum
    cdfTable <- t(apply(normTable, 1, cumsum))
    
    #wasserstein distance
    wd <- dist(cdfTable, method = "manhattan")
    
    ## summarize distances ##
    
    #table of all unique pairs of profiles
    d <- data.table(Protein_1 = ids)
    d[, `:=`(id1 = 1L, id2 = .I)]
    setkey(d, id1, id2)
    distances <- data.table::foverlaps(d, d, type = "within", which = TRUE)[xid != yid]
    distances <- data.table::setDT(list(d$Protein_1[distances$xid], d$Protein_1[distances$yid]))
    
    #get Wasserstein distance for each pairing
    distances$distance <- as.numeric(wd)
    
    #rename columns
    colnames(distances) <- c("id1", "id2", "distance")
    
  }
  
  #euclidean distance
  if(metric == "euclidean"){
    
    #normalize profiles to sum to 1
    normTable <- t(apply(cofracTable, 1, function(profile){
      return(profile/sum(profile))
    }))
    
    #euclidean distance
    eucl <- dist(normTable, method = "euclidean")
    
    #summarize distances#
    
    #table of all unique pairs of profiles
    d <- data.table(Protein_1 = ids)
    d[, `:=`(id1 = 1L, id2 = .I)]
    setkey(d, id1, id2)
    distances <- data.table::foverlaps(d, d, type = "within", which = TRUE)[xid != yid]
    distances <- data.table::setDT(list(d$Protein_1[distances$xid], d$Protein_1[distances$yid]))
    
    #get Wasserstein distance for each pairing
    distances$distance <- as.numeric(eucl)
    
    #rename columns
    colnames(distances) <- c("id1", "id2", "distance")
  }
  
  #pearson correlation
  if(metric == "pearson"){
    #transpose matrix
    cofracTable <- t(cofracTable)
    colnames(cofracTable) <- (ids)
    
    #compute all pairwise pearson correlations
    cors <- cor(cofracTable, method = "pearson")
    
    #summarize distances
    distances <- data.frame(id1=rownames(cors)[row(cors)], id2=colnames(cors)[col(cors)], distance=c(cors))
    
    #rename columns
    colnames(distances) <- c("id1", "id2", "distance")
    
    #convert pearson similarity to distance measure
    distances$distance <- 1 - distances$distance
  }
  
  #return matrix with pairwise distances
  return(as.data.frame(distances))
}
