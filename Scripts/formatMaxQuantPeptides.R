#' @name formatMaxQuantPeptides
#' @title Format MaxQuant peptides file
#'
#' @description Formats TMT-multiplexed peptides.txt file from MaxQuant output to
#'  programming-friendly peptide-by-fraction tables.
#'
#' @param file_path Character vector giving the list of local paths
#'  to the peptides.txt files. One file per replicate.
#' @param conditions Character vector with names of experimental conditions.
#'
#' @return List of lists. One list per replicate named \code{replicate_X}.
#' Each replicate list consists of one peptide-by-fraction data frame per condition.
#' The first column of data frame has the peptide sequence,
#' the second column has the gene name,
#' and the following columns contain the corrected reporter intensities for the fractions.
#' Only the peptides that were detected in all replicates are retained.
#'
#' @examples
#' peptides <- formatMaxQuantPeptides(c("Documents/peptides_rep1.txt", "Documents/peptides_rep2.txt"),
#'  conditions = c("Brain", "Lung", "Liver"))
#'
#' @import dplyr data.table
NULL
#' @export
formatMaxQuantPeptides <- function(file_path, conditions){
  #iterate over files (one per replicate)
  peps <- lapply(file_path, function(file){
    #read in MaxQuant file
    peps <- as.data.frame(data.table::fread(file))

    #make sure correct columns are present
    if(sum(c("Sequence", "Proteins") %in% colnames(peps)) < 2){
      stop("Incorrect format. Please input peptides.txt file from MaxQuant output.")
    }

    #remove contaminants
    peps <- peps[-grep("CON", peps$Proteins),]

    #remove reverse decoys in peptides
    peps <- peps[-grep("REV", peps$Sequence),]

    #keep relevant columns only (peptide/protein name, reporter corrected intensity values)
    peps <- peps[,c("Sequence", "Gene names",
                    grep(pattern = "corrected", x = colnames(peps), value = T))][,-c(3:(length(conditions)+2))]

    #create empty list to store data for each condition
    conds <- vector(mode = "list", length = length(conditions))
    names(conds) <- conditions

    #create peptide/protein table list for easier downstream processing
    #loop over conditions
    for(i in seq_along(conditions)){
      #get relevant condition
      tmp <- peps[ , c(1, 2, grep(pattern = paste0("Reporter intensity corrected ",i-1,"[[:blank:]][[:alpha:]]"), x = colnames(peps)))]

      #rename column names
      colnames(tmp) <- c("sequence", "protein",
                         gsub(pattern = paste0("Reporter intensity corrected ", i-1, "[[:blank:]][[:alpha:]]"), replacement = "", x = colnames(tmp)[-c(1,2)]))

      #reorder columns by fractions
      tmp <- tmp[,c(1, 2, order(as.numeric(colnames(tmp)[-c(1,2)]))+2)]

      #set rownames to unique sequences
      rownames(tmp) <- tmp[,1]

      #add to list
      conds[[i]] <- tmp
    }

    #change variable names and remove temporary variables
    peps <- conds
    rm(conds); rm(tmp)

    #remove nameless proteins
    peps <- lapply(peps, function(curCond){return(curCond[curCond$protein != "",])})
  })

  #name lists according to corresponding replicate
  names(peps) <- paste0("replicate_", seq(length(peps)))

  #if multiple replicates present
  if(length(peps) > 1){
    ## retain overlapping peptides only ##
    #get all peptides detected in each replicate
    allPeps <- lapply(peps, function(curRep){
      return(unlist(lapply(curRep, "[", 1)))
    })

    #find overlapping peptides
    commonPeps <- Reduce(intersect, allPeps)

    #subset tables to overlapping peptides only
    peps <- lapply(peps, function(curRep){
      return(lapply(curRep, function(curCond){
        return(curCond[commonPeps,])
      }))
    })
  }
  return(peps)
}
