## Merge peptide co-fractionation profiles ##

#required packages
library(dplyr) #data wrangling
library(reshape2) #data wrangling

# input: Profile-by-fraction matrix, first column has peptide IDs, second column has protein IDs and must be named 'protein'.
# input: Method for merging peptide profiles.
# output: One protein-by-fraction matrix.
mergePeptideProfiles <- function(pepTable, method = c("sum", "average")){
  #convert data to long format
  pepTable <- reshape2::melt(pepTable)
 
  #user message
  print("Creating protein profiles from peptide profiles...")
  
  #sum up peptide profiles belonging to each protein
  if(method == "sum"){
    #sum peptides in each fraction
    proteinTable <- pepTable %>%
      group_by(protein, variable) %>%
      summarise(intensity = sum(value))
    
    #convert to wide format
    proteinTable <- reshape2::dcast(proteinTable, protein ~ variable)
    
    #remove proteins with zero intensity values
    proteinTable <- proteinTable[rowSums(proteinTable[,-1]) > 0, ]
    
    return(proteinTable)
  }
  
  #average peptide profiles belonging to each protein
  if(method == "average"){
    #average peptides in each fraction
    proteinTable <- pepTable %>%
      group_by(protein, variable) %>%
      summarise(intensity = mean(value))
    
    #convert to wide format
    proteinTable <- reshape2::dcast(proteinTable, protein ~ variable)
    
    #remove proteins with zero intensity values
    proteinTable <- proteinTable[rowSums(proteinTable[,-1]) > 0, ]
    
    return(proteinTable)
  }
}
