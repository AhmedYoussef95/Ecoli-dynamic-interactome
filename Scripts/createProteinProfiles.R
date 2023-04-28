#' @name createProteinProfiles
#' @title Create protein profiles from peptide-level CF/MS data
#'
#' @description Collapse peptide-level dataset to the corresponding protein-level intensities.
#'
#' There are three steps involved:
#' * ***Filter outlier peptides***: Each peptide is mapped to a particular
#'  protein by the peptide mapping software like *MaxQuant*
#'  in an upstream step. It is expected that the peptides that map to the same
#'  protein should have similar profiles to each other, and as such any peptide
#'  that deviates significantly from its group is likely a faulty measurement.
#'  Two strategies for identifying these measurements are described below.
#'
#'  * ***Select top-intensity peptides***:  A common strategy for creating
#'  protein profiles from peptide profiles is to average or sum the intensities
#'   of the top two or three high-intensity peptides to represent the final
#'    protein profile (e.g. Silva et. al, 2005). This function gives the user
#'    the option to perfrm such filtering based on any given number of peptides.
#'
#'  * ***Collapse peptide intensities to protein-level intensities***:
#'  Sum or average sibling peptide intensities across fractions to create
#'  the corresponding protein-level profiles.
#'
#' @param peps List of lists. One list per replicate.
#'  Each replicate list consists of one peptide/protein-by-fraction data frame
#'  per condition.
#' @param cores Integer indicating the number of computer cores to use for
#' parallel computation of distances.
#' @param filterOutliers Character variable indicating type of sibling peptide
#'  outlier filtering to be applied. The three options are:
#'  * "clustering": Perform average-linkage hierarchical clustering on the
#'   sibling peptides based on their similarity to each other, split the
#'    resulting dendrogram into two clusters, and retain the peptides belonging
#'     to the larger cluster as being representative of the protein.
#'  * "quantile": Compute each peptide’s average similarity to its
#'   ‘sibling peptides’ and filter out the peptides beyond a quantile-based
#'    threshold, which by default is 0.95. This means that the 95% most similar
#'     peptides would be retained for downstream analysis.
#'  * "none" *or any other value*: No peptide outlier filtering.
#'   Any value aside from the above two will not apply any filtering
#'    to the data.
#'    Default method is "clustering".
#' @param threshold Numeric between 0 and 1 indicating threshold to filter
#'  outlier peptides based upon. Only relevant if quantile filtering is
#'  chosen. Default is 0.95.
#' @param topN Integer indicating how many peptides to retain per protein.
#'  Selected peptides are those with the highest cross-fraction summed
#'  intensities. Default is NA which retains all peptides.
#' @param distanceMetric Character indicating which distance metric to use to
#' compute sibling peptide similarity based upon. Choices are Wasserstein
#'  distance, Euclidean distance, and Pearson distance (1 - Pearson R2).
#'  Default is Pearson.
#' @param method Character specifying how to merge peptide intensities to
#'   protein intensities. One of "sum" or "average". Default is "average".
#'
#' @return List of lists. One list per replicate named \code{replicate_X}.
#' Each replicate list consists of one protein-by-fraction data frame per condition.
#' The first column of data frame has the protein ID,
#' and the following columns contain the intensities for each fraction.
#' Only the peptides that were detected in all replicates are retained.
#'
#' @examples
#' proteins <- createProteinProfiles(peps, topN = 3)
#'
#' @import dplyr data.table data.table reshape2 igraph
NULL
#' @export
#'
createProteinProfiles <- function(peps,
                                  cores = 1,
                                  filterOutliers = c("clustering", "quantile", "none"),
                                  threshold = 0.95,
                                  topN = NA,
                                  distanceMetric = c( 'pearson', 'wasserstein', 'euclidean'),
                                  method = c("average", "sum")){
  #filter outlier peptides
  if(filterOutliers %in% c("quantile", "clustering")){
    peps <- lapply(peps, function(curReplicate){
      return(lapply(curReplicate, filterOutlierPeptides,
                    metric = distanceMetric, method = filterOutliers,
                    cores = cores, threshold = NULL))
    })
  }

  #select top N intensity proteins per protein
  if(is.numeric(topN)){
    peps <- lapply(peps, function(curReplicate){
      return(lapply(curReplicate, selectTopPeptides, n = topN))
    })
  }

  #merge peptide profiles into protein profiles
  proteins <- lapply(peps, function(curReplicate){
    return(suppressMessages(
      lapply(curReplicate, mergePeptideProfiles, method = method))) })

  return(proteins)
}

## Helper function: compute distances between sibling peptides ##
# input: profile-by-fraction matrix, first column has peptide IDs,
# second column has protein IDs and must be named 'protein'
# output: list of sibling peptide distance matrices, one per protein
siblingDistances <- function(pepTable, metric = c('wasserstein', 'euclidean', 'pearson'), cores = 6){

  #split input table into tables of sibling peptides
  profiles <- split(pepTable, with(pepTable, protein), drop=TRUE)

  #compute all pairwise distances for each group of sibling peptides
  sibDistances <- pblapply(profiles, function(siblings){
    #if single-peptide protein return a distance of 0
    if(nrow(siblings) == 1)
      sibs <- data.frame(id1 = siblings[,1],  id2 = siblings[,1], distance = 0)
    else
      sibs <- computeDistances(cofracTable = siblings, metric = metric)
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
                                  metric = c('pearson', 'euclidean', 'wasserstein'),
                                  method = c("clustering", "quantile"),
                                  threshold = 0.95,
                                  cores){
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
    #print("Computing average agreement of each peptide with its siblings...")
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
    #print("Filtering outlier peptides using clustering...")
    #iterate over proteins and decide which peptides to retain
    clusteredPeps <- pblapply(sibDistances, function(siblings){
      #remove any NA values (non-reproducible peptides)
      siblings <- siblings[complete.cases(siblings),]

      #if single-peptide protein return the peptide
      if(length(union(siblings$id1, siblings$id2)) <= 1){
        return(siblings$id1)
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
  pepSums <- data.table::data.table(cbind(pepTable[,c(1,2)], pepSums), key = "protein")

  #in case of replicates, keep max sum across replicates for each peptide
  pepSums <- pepSums %>%
    dplyr::group_by_at(1) %>%
    dplyr::filter(pepSums == max(pepSums))

  #find top N peptides for each protein
  topPeps <- pepSums %>%
    dplyr::arrange(desc(pepSums)) %>%
    dplyr::group_by(protein) %>%
    dplyr::slice(1:n)

  #retain only the top N peptides for each protein
  filteredPeps <- pepTable[unlist(pepTable[,1]) %in% unlist(topPeps[,1]),]

  return(filteredPeps)
}

# input: Profile-by-fraction matrix, first column has peptide IDs, second column has protein IDs and must be named 'protein'.
# input: Method for merging peptide profiles.
# output: One protein-by-fraction matrix.
mergePeptideProfiles <- function(pepTable, method = c("average", "sum")){
  #convert data to long format
  suppressMessages(pepTable <- reshape2::melt(pepTable))

  #sum up peptide profiles belonging to each protein
  if(method == "sum"){
    #sum peptides in each fraction
    proteinTable <- pepTable %>%
      dplyr::group_by(protein, variable) %>%
      dplyr::summarise(intensity = sum(value))

    #convert to wide format
    proteinTable <- reshape2::dcast(proteinTable, protein ~ variable)

    #remove proteins with zero intensity values
    proteinTable <- proteinTable[rowSums(proteinTable[,-1]) > 0, ]

    #set rownames to protein IDs
    rownames(proteinTable) <- proteinTable$protein

    return(proteinTable)
  }

  #average peptide profiles belonging to each protein
  if(method == "average"){
    #average peptides in each fraction
    proteinTable <- pepTable %>%
      dplyr::group_by(protein, variable) %>%
      dplyr::summarise(intensity = mean(value))

    #convert to wide format
    proteinTable <- reshape2::dcast(proteinTable, protein ~ variable)

    #remove proteins with zero intensity values
    proteinTable <- proteinTable[rowSums(proteinTable[,-1]) > 0, ]

    #set rownames to protein IDs
    rownames(proteinTable) <- proteinTable$protein

    return(proteinTable)
  }
}
