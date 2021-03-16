# Script for pre-processing MS co-frac peptide data #

#required packages
library(forecast) #smoothing data with moving average
library(pbapply) #parallelize operations

## Normalize data ##

#function to normalize column-wise
normFrac <- function(pepTable){
  return(cbind(pepTable[,1:2], t(t(pepTable[,3:98])/colSums(pepTable[,3:98]))))
}

## Smooth signals with moving average ##

#function for smoothing all profiles in a given table
smoothProfile <- function(pepTable, windowSize = 4){
  #iterate over all rows of the table
  smoothed <- apply(pepTable[,-c(1,2)], 1, function(profile){
    toSmooth <- c(rep(0,round(windowSize/2)), t(profile), rep(0,round(windowSize/2)))
    return(forecast::ma(toSmooth, windowSize)[3:98])
  })
  
  return(cbind(pepTable[,c(1,2)], t(smoothed)))
}
