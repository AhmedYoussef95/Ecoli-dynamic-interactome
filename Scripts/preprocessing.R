#' @name preprocessing
#' @title Preprocess CF/MS profiles
#'
#' @description  Log-transform and smooth CF/MS elution profiles.
#'
#' @param profiles List of lists. One list per replicate.
#'  Each replicate list consists of one peptide/protein-by-fraction data frame
#'  per condition.
#' @param logTransform Character variable indicating type of log-transformation
#'  to be applied. The four options are:
#'  * "log2": Log2-transform profiles. A pseudocount of 1 is added to each
#'   value prior to log-transformation to avoid computing infinity values.
#'  * "CLR": Centered log-transform. This is achieved by subtracting the mean
#'   value for each peptide/protein in each fraction (i.e. across conditions)
#'    after log2-transformation.
#'  * "normCLR": Normalized centered log-transform. Approach adopted from Drew
#'   et al. (2017) which normalizes the data within each fraction, across all
#'    experimental conditions but within each replicate separately,
#'    followed by a centered log2-transformation.
#'  * "none" *or any other value*: No log-transformation.
#'   Any value aside from the above three will not apply any transformation
#'   to the data.
#' @param smooth Logical variable indicating whether to smooth CF/MS profiles
#'  using moving average.
#'  Default is \code{FALSE}.
#' @param smoothingFractions Integer indicating the number of fractions
#'  to be considered for smoothing the profiles using a moving average window.
#'  Default is 4 fractions. Irrelevant if \code{smooth = TRUE}.
#'
#' @return List of lists. Same structure as the \code{profiles} input list.
#'
#' @examples
#' peptides <- preprocessing(profiles = profiles, logTransform = "log",
#'  smooth = TRUE, smoothingFractions = 4)
#'
#' @import dplyr data.table reshape2 forecast
NULL
#' @export
preprocessing <- function(profiles,
                          logTransform = c('log2','CLR', 'normCLR','none'),
                          smooth = c(FALSE, TRUE),
                          smoothingFractions = 4){
  #log-transform
  if(logTransform %in% c("log2", "CLR", "normCLR")){
    if(logTransform == "normCLR"){
      #normalize profiles
      profiles <- lapply(profiles, function(curReplicate){
        #label conditions
        curReplicate <- lapply(names(curReplicate), function(curCond){
          curReplicate[[curCond]]$cond <- curCond
          return(curReplicate[[curCond]])
        })

        #normalize across each fraction
        suppressMessages(normprofiles <- normFrac(curReplicate))

        return(normprofiles)
      })
    }

    #log2-transform
    profiles <- lapply(profiles,function(curRep){
      return(lapply(curRep, function(curCond){
        #store profile (protein/peptide) IDs in variable
        ids <- dplyr::select_if(curCond, is.character)

        #keep intensity values only
        cofracTable <- dplyr::select_if(curCond, is.numeric)

        #log2-transformation with pseudocount of 1
        cofracTable <- log2(cofracTable + 1)

        return(cbind(ids, cofracTable))
      }))
    })

    #centered log-transform
    if(logTransform %in% c("CLR", "normCLR")){
      profiles <- lapply(profiles, function(curReplicate){
        return(lapply(curReplicate, function(curCond){
          #store profile (protein/peptide) IDs in variable
          ids <- dplyr::select_if(curCond, is.character)

          #keep intensity values only
          cofracTable <- dplyr::select_if(curCond, is.numeric)

          #subtract mean in each fraction
          cofracTable <- apply(cofracTable, 2, function(x) return(x - mean(x)))
          return(cbind(ids, cofracTable))
        }))
      })
    }
  }

  #smooth profiles
  if(smooth){
    profiles <- lapply(profiles, function(curReplicate){
      return(lapply(curReplicate, smoothProfile, windowSize = smoothingFractions))
    })
  }

  return(profiles)
}

#' Normalize CF/MS profiles within fractions
#'
#' Helper function to normalize peptide/protein CF/MS profiles within each
#' fraction and across experimental conditions.
#'
#' @param condTables List of dataframes. One dataframe per experimental
#'  condition. The first columns have the protein ID, and peptide sequence
#'  if the input is peptide-level data, and the subsequent columns have the
#'  per-fraction intensities.
#'
#' @return List of dataframes after normalization. Same structure as the
#'  \code{condTables} input list.
normFrac <- function(condTables){
  #merge conditions into one table
  condTables <- do.call("rbind", condTables)

  #convert to long format
  condTables <- reshape2::melt(condTables)

  #normalize intensities of each peptide in each fraction across conditions
  condTables <- condTables %>%
    dplyr::group_by(variable) %>%
    dplyr::mutate(value = value / (sum(value) + 0.1))

  #separate peptides by condition
  condTables <- split(condTables, condTables$cond)

  #re-convert to wide format
  if("sequence" %in% colnames(condTables))
    condTables <- lapply(condTables, function(x){
      return(reshape2::dcast(x, sequence + protein ~ variable))})
  else
    condTables <- lapply(condTables, function(x){
      return(reshape2::dcast(x, protein ~ variable))})

  return(condTables)
}

#' Smooth signals with moving average
#'
#' Helper function to smooth peptide/protein CF/MS profiles using a
#' moving average approach. The *forecast* package is used in this function.
#'
#' @param condTables List of dataframes. One dataframe per experimental
#'  condition. The first columns have the protein ID, and peptide sequence,
#'  if the input is peptide-level data, and the subsequent columns have the
#'  per-fraction intensities.
#' @param windowSize Integer indicating the number of fractions
#'  to be considered for smoothing the profiles using a moving average window.
#'
#' @return List of dataframes after smoothing. Same structure as the
#'  \code{condTables} input list.
#function for smoothing all profiles in a given table
smoothProfile <- function(condTables, windowSize){
  #store profile (protein/peptide) IDs in variable
  ids <- dplyr::select_if(condTables, is.character)

  #keep intensity values only
  cofracTable <- dplyr::select_if(condTables, is.numeric)

  #iterate over all rows of the table
  smoothed <- apply(cofracTable, 1, function(profile){
    toSmooth <- c(rep(0,round(windowSize/2)), t(profile), rep(0,round(windowSize/2)))
    return(forecast::ma(toSmooth, windowSize)[3:98])
  })

  return(cbind(ids, t(smoothed)))
}
