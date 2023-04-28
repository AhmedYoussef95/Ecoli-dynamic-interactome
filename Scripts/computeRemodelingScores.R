#' @name computeRemodelingScores
#' @title Compute PPI remodeling scores
#'
#' @description Compute conditional change in PPI similarity scores relative
#'  to reference condition.
#'  
#' @param conditionalSimilarities Data frame with PPI conditional similarity
#'  scores. Output of \code{computeConditionalSimilarities}.
#' 
#' @param referenceCondition Either name of reference condition or index of
#' reference condition column in \code{conditionalSimilarities}
#' 
#' @return Data frame. Remodeling score for each PPI in each condition.
#' The higher the score, the more likley the PPI is disrupted.
#' 
#' @examples
#' remodelingScores <- computeRemodelingScores(conditionalSimilarities,
#'  referenceCondition = "LB")
#'
#' @import dplyr magrittr WGCNA reshape2
NULL
#' @export
computeRemodelingScores <- function(conditionalSimilarities,
                                           referenceCondition,){
  
  #compute remodeling scores as difference between each condition and reference condition
  ppiRemodelScores <- lapply(colnames(conditionalSimilarities)[-c(1,2)],
                             function(curCond){
    return(conditionalSimilarities[,referenceCondition] -
             conditionalSimilarities[,curCond]) 
  }) %>% 
    setNames(conditions) %>% 
    cbind(conditionalSimilarities[,c(1,2)], .) %>% 
    select(-referenceCondition)

  return(ppiRemodelScores)
}