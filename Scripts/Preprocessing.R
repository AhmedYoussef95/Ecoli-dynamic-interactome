# Script for pre-processing MS co-frac peptide data #

#required packages
library(data.table) #fast reading of input data
library(pbapply) #parellelize operations
library(dplyr) #data wrangling
library(tidyr) #data wrangling
library(forecast) #smoothing data with moving average

## 1. Read in and format data ##

#set working directory to folder containing data
setwd('./Raw data/')

#read in peptide-level data for both replicates
peps1 <- as.data.frame(fread('peptides.txt'))
peps2 <- as.data.frame(fread('peptides_rep2.txt'))

#remove contaminants
peps1 <- peps1[-grep("CON", peps1$Proteins),]
peps2 <- peps2[-grep("CON", peps2$Proteins),]

#find the common proteins between the 2 replicates
commonProts <- intersect(peps1$`Gene names`, peps2$`Gene names`)

#keep only the peptides mapping to common proteins
peps1 <- peps1[peps1$`Gene names` %in% commonProts,]
peps2 <- peps2[peps2$`Gene names` %in% commonProts,]

#find the common peptides between the 2 replicates
commonPeps <- intersect(peps1$Sequence, peps2$Sequence)

#keep the common peptides only (29,630 peptides)
peps1 <- peps1[peps1$Sequence %in% commonPeps,]
peps2 <- peps2[peps2$Sequence %in% commonPeps,]

#keep relevant columns only (peptide/protein name, reporter corrected intenisty values)
peps1 <- peps1[,c(1, 39, grep(pattern = "corrected", x = colnames(peps1)))[-c(3:13)]]
peps2 <- peps2[,c(1, 39, grep(pattern = "corrected", x = colnames(peps2)))[-c(3:13)]]

#vector with names of experimental conditions
conditions <- c("LB","glucose","galactose","xylose","aa","42c","sp","nzg","max","ana","combined")

#create empty list to store data for each condition (2 replicates separate)
conds1 <- conds2 <- vector(mode = "list", length = 11)
names(conds1) <- names(conds2) <- conditions

#create peptide/protein table list for easier downstream processing
#loop over conditions
for(i in seq_along(conditions)){
  #get relevant condition
  tmp1 <- peps1[, c(1,2,grep(pattern = paste0("Reporter intensity corrected ",i-1," F"), x = colnames(peps1)))]
  tmp2 <- peps2[, c(1,2,grep(pattern = paste0("Reporter intensity corrected ",i-1," F"), x = colnames(peps2)))]
  
  #rename fraction names
  colnames(tmp1) <- colnames(tmp2) <- c("sequence","gene",gsub(paste0("Reporter intensity corrected ",i-1," F"),"", colnames(tmp1)[3:98]))
  
  #reorder columns by fractions
  tmp1 <- tmp1[,c(1,2,order(as.numeric(colnames(tmp1)[3:98]))+2)]
  tmp2 <- tmp2[,c(1,2,order(as.numeric(colnames(tmp2)[3:98]))+2)]
  
  #add to list
  conds1[[i]] <- tmp1
  conds2[[i]] <- tmp2
}

#change variable names and remove temporary variables
peps1 <- conds1; peps2 <- conds2
rm(conds1); rm(conds2); rm(tmp1); rm(tmp2)

## 2. Normalize data ##

#function to normalize
normFrac <- function(pepTable){return(cbind(pepTable[,1:2], t(t(pepTable[,3:98])/colSums(pepTable[,3:98]))))}

#normalize both replicates
normPeps1 <- lapply(peps1, normFrac)
normPeps2 <- lapply(peps2, normFrac)

## 3. Smooth signals with moving average ##

#function for smoothing all profiles in a given table
smoothProfile <- function(table){
  #iterate over all rows of the table
  smoothed <- apply(table[,-c(1,2)], 1, function(profile){
    noise <- runif(n = 2, min = 0, max = 1e-10)
    toSmooth <- c(rep(0,2), t(profile), rep(0,2))
    return(forecast::ma(toSmooth, 4)[3:98])
  })
  
  return(cbind(table[,c(1,2)], t(smoothed)))
}

#smooth all peptide profiles (~10s each)
smoothedPeps1 <- pblapply(normPeps1, smoothProfile, cl = 6)
smoothedPeps2 <- pblapply(normPeps2, smoothProfile, cl = 6)
