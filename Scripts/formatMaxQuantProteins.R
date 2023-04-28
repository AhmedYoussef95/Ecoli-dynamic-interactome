#' @name formatMaxQuantProteins
#' @title Format MaxQuant proteinGroups.txt file
#'
#' @description Formats TMT-multiplexed proteinGroups.txt file from MaxQuant output to
#'  programming-friendly protein-by-fraction tables.
#'
#' @param file_path Character vector giving the list of local paths
#'  to the proteinGroups.txt files. One file per replicate.
#' @param conditions Character vector with names of experimental conditions.
#'
#' @return List of lists. One list per replicate named \code{replicate_X}.
#' Each replicate list consists of one protein-by-fraction data frame per condition.
#' The first column of data frame has the protein ID as stored in the 'Gene names' column,
#' and the following columns contain the corrected reporter intensities for the fractions.
#' Only the proteins that were detected in all replicates are retained.
#'
#' @examples
#' proteins <- formatMaxQuantProteins(c("Documents/proteinGroups_rep1.txt", "Documents/proteinGroups_rep2.txt"),
#'  conditions = c("Brain", "Lung", "Liver"))
#'
#' @import dplyr data.table
NULL
#' @export
formatMaxQuantProteins <- function(file_path, conditions){
  #iterate over files (one per replicate)
  prots <- lapply(file_path, function(file){
    #read in MaxQuant file
    prots <- as.data.frame(data.table::fread(file))

    #make sure correct columns are present
    if("Sequence" %in% colnames(prots)){
      stop("Incorrect format. Please input proteinGroups.txt file from MaxQuant output.")
    }

    #remove contaminants
    prots <- prots[-grep("^CON", prots$`Protein IDs`),]

    #keep relevant columns only (peptide/protein name, reporter corrected intensity values)
    tmp <- prots[,c("Gene names",
                      grep(pattern = "corrected", x = colnames(prots), value = T))][,-c(3:(length(conditions)+2))]

    #create empty list to store data for each condition
    conds <- vector(mode = "list", length = length(conditions))
    names(conds) <- conditions

    #create peptide/protein table list for easier downstream processing
    #loop over conditions
    for(i in seq_along(conditions)){
      #get relevant condition
      tmp <- prots[ , c(1, grep(pattern = paste0("Reporter intensity corrected ",i-1,"[[:blank:]][[:alpha:]]"), x = colnames(prots)))]

      #rename column names
      colnames(tmp) <- c("protein",
                         gsub(pattern = paste0("Reporter intensity corrected ", i-1, "[[:blank:]][[:alpha:]]"), replacement = "", x = colnames(tmp)[-1]))

      #reorder columns by fractions
      tmp <- tmp[,c(1, order(as.numeric(colnames(tmp)[-1]))+1)]

      #set rownames to unique IDs
      rownames(tmp) <- tmp[,1]

      #add to list
      conds[[i]] <- tmp
    }

    #change variable names and remove temporary variables
    prots <- conds
    rm(conds); rm(tmp)

    #remove nameless proteins
    prots <- lapply(prots, function(curCond){return(curCond[curCond$protein != "",])})
  })

    #name lists according to corresponding replicate
    names(prots) <- paste0("replicate_", seq(length(prots)))

    #if multiple replicates present
    if(length(prots) > 1){
      ## retain overlapping proteins only ##
      #get all proteins detected in each replicate
      allprots <- lapply(prots, function(curRep){
        return(unlist(lapply(curRep, "[", 1)))
      })

      #find overlapping proteins
      commonprots <- Reduce(intersect, allprots)

      #subset tables to overlapping proteins only
      prots <- lapply(prots, function(curRep){
        return(lapply(curRep, function(curCond){
          return(curCond[commonprots,])
        }))
      })
    }
    return(prots)
}
