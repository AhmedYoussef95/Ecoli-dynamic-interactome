## Script for comparing PPIs to external refeerence databases ##
#required packages
library(readxl) #working with excel files
library(PRROC) #ROC curve

#read in and format 4 reference Ecoli PPIs
assembleRefNetwork <- function(ppi){
  #list of proteins in our dataset
  ourProts <- union(ppi[,1], ppi[,2])
  
  ## Y2H binary PPIs
  
  #read in excel file
  y2h <- as.data.frame(readxl::read_excel("/Users/Ahmed/Documents/Emili_Lab/Aim 1/Raw data/binary_PPIs_rajogapala_2014.xls",
                                          skip = 2, col_names = c("A","B","id1","id2","source")))
  
  #keep only experimentally-validated interactions
  y2h <- y2h[y2h$source!="LIT", 1:2]
  
  #keep gene names only
  y2h[,1] <- gsub("_.*","", y2h[,1]); y2h[,2] <- gsub("_.*","", y2h[,2])
  
  #keep only the pairs where both proteins are in our co-frac dataset
  y2h <- y2h[y2h[,1]%in%ourProts & y2h[,2]%in%ourProts,]
  
  ## Hu et al. AP/MS Complexes
  
  #read in reference complexes
  APMS <- readxl::read_excel('/Users/Ahmed/Documents/Emili_Lab/Aim 1/Raw data/Ecoli_reference_complexes.xlsx', sheet = 2, skip = 1)[,2:3]
  colnames(APMS) <- c("A","B")
  
  #store complexes as list
  APMS <- lapply(unlist(APMS[,2]), function(x) unlist(strsplit(x, split = ", ")))
  
  #for each complex get all possible member pairs and combine into one big table
  APMS <- do.call("rbind",lapply(APMS, function(x) t(combn(x,2))))
  
  #keep only the pairs where both proteins are in our co-frac dataset
  APMS <- APMS[APMS[,1]%in%ourProts & APMS[,2]%in%ourProts,]
  
  #re-adjust column names
  colnames(APMS) <- c("A","B")
  
  ## STRING
  
  #read in binary PPIs and supporting info
  string <- read.delim("/Users/Ahmed/Documents/Emili_Lab/Aim 1/Raw data/string_interactions.txt", sep = " ", stringsAsFactors = F)
  stringAnno <- read.delim("/Users/Ahmed/Documents/Emili_Lab/Aim 1/Raw data/string_protein_info.txt", stringsAsFactors = F)
  
  #filter annotation table
  stringAnno <- stringAnno[stringAnno$preferred_name %in% ourProts,]
  
  #keep only the pairs where both proteins are in our co-frac dataset
  string <- string[string[,1]%in%stringAnno$protein_external_id & string[,2]%in%stringAnno$protein_external_id,]
  
  #keep high-confidence pairs only
  string <- string[string$combined_score >= 700, 1:2]
  
  #replace gene IDs in PPI table with gene names
  string <- t(pbapply(string, 1, function(x){
    x[1] <- stringAnno$preferred_name[stringAnno$protein_external_id==x[1]]
    x[2] <- stringAnno$preferred_name[stringAnno$protein_external_id==x[2]]
    return(x)
  }, cl = 6))
  
  #column names
  colnames(string) <- c("A","B")
  
  ## BioGRID
  
  #read in data
  biogrid <- read.delim("/Users/Ahmed/Documents/Emili_Lab/PPI Databases/BioGRID/BIOGRID-ORGANISM-Escherichia_coli_K12_W3110-3.5.186.tab3.txt")[,8:9]
  
  #keep only the pairs where both proteins are in our co-frac dataset
  biogrid <- setNames(biogrid[biogrid[,1] %in% ourProts & biogrid[,2] %in% ourProts,], c("A", "B"))
  
  ## EPIC-predicted network
  epic <- as.data.frame(fread("/Users/Ahmed/Documents/Emili_Lab/Aim 1/EPIC results/pooled/results/Out.pred.txt",
                              header = FALSE, stringsAsFactors = FALSE))[,1:2]
  
  #convert UniProt IDs to gene names
  uniprot <- read.csv("/Users/Ahmed/Documents/Emili_Lab/Aim 1/Processed data/uniprot_mapping.csv",
                      header = FALSE, stringsAsFactors = FALSE, skip = 1)
  uniprot <- setNames(uniprot[,1], uniprot[,2])
  epic[,1] <- uniprot[unlist(epic[,1])]; epic[,2] <- uniprot[unlist(epic[,2])]
  
  #column names
  colnames(epic) <- c("A","B")
  
  #final list of networks
  refNetworks <- list(y2h = y2h, APMS = APMS, string = string, biogrid = biogrid, epic = epic)
  
  #pool into one network
  refPPI <- do.call("rbind", refNetworks)
  
  #sort network alphabetically to maintain same order
  refPPI <- t(apply(refPPI, 1, sort))
  
  #remove duplicate interactions
  refPPI <- unique(refPPI)
  
  #remove any self-interactions
  refPPI <- as.data.frame(refPPI[refPPI[,1] != refPPI[,2],])
  
  #create key for each PPI
  ppi$key <- unlist(apply(ppi[,1:2], 1, paste, collapse = " "))
  refPPI$key <- unlist(apply(refPPI, 1, paste, collapse = " "))
  
  #get ppi scores for this reference PPI
  refPPI <- ppi[ppi$key %in% refPPI$key,]
  
  #retain only PPIs that exhibit greater-than-average similarity in at least one condition in our dataset 
  cutoff <- quantile(as.matrix(ppi[,-c(1,2, ncol(ppi))]), 0.01)
  refPPI <- refPPI[unlist(apply(refPPI[,-c(1,2,ncol(refPPI))],1 , min)) < cutoff,]
  
  return(refPPI[,-ncol(refPPI)])
}
