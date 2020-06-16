#' Convert a data frame into an object of class \code{tfa_metadata}
#'
#' \code{as.tfa_metadata} checks the contents of a metadata file for the presence of group and instance IDs,
#' and it convert the data frame as an object of class \code{tfa_metadata}.
#'
#' @param dataframe a dataframe containing metadata for a factor analysis of ancient DNA. The meta
#' data could contain information on individual ID's, group ID's, country of origin, sample age, etc. Columns
#' with names 'Instance.ID' and 'Group.ID' are mandatory.
#' @return an object of class \code{tfa_metadata}
#' @export
#' @examples
#' library(tfa)
#'
#' # Ancient DNA from Bronze Age Great Britain samples
#' # including Steppe (Yamnaya), hunter gatherers, and early farmers from Anatolia
#' data(england_ba)
#' attach(England_BA)
#' metadata <- as.tfa_metadata(meta)
#' detach(England_BA)
#' @references François, O., Jay, F. (2020). Factor analysis of ancient DNA samples.
#' @seealso \code{\link{England_BA}}, \code{\link{ancestry_coefficients}}
as.tfa_metadata <- function(dataframe){

  if (class(dataframe) != "data.frame"){
   stop("Object dataframe not of class data.frame")
  }
  cnames <- colnames(dataframe)
  if ("Group.ID" %in% cnames == FALSE){
    stop("dataframe has no column named Group.ID")
  }
  if (!is.character(dataframe[,"Group.ID"])){
    stop("Group.ID is not character")
  }
  if ("Instance.ID" %in% cnames == FALSE){
    stop("dataframe has no column named Instance.ID")
  }
  obj <- dataframe
  class(obj) <- "tfa_metadata"
  return(obj)
}

#' Compute ancestry coefficients from specified source populations
#'
#' \code{ancestry_coefficients} compute the proportions of shared ancestry from specified source populations
#'
#' @param model an object of class \code{tfa} with the same number of individuals as in metadata.
#' @param metadata an object of class \code{tfa_metadata} containing metadata on Group.IDs.
#' @param source a vector of ancestral group ID's from the list of \code{Group.IDs} (character).
#' @param target a vector of target group ID's from the list of \code{Group.IDs} (character).
#' @param individual a logical indicating whether individual ancestry coefficients
#' should be computed.
#' @return a matrix with coefficients equal to the ancestry proportions of each target population
#' from each source.
#' @export
#' @examples
#' library(tfa)
#'
#' # Ancient DNA from Bronze Age Great Britain samples
#' # including Steppe (Yamnaya), hunter gatherers, and early farmers from Anatolia
#' data(england_ba)
#' attach(England_BA)
#' metadata <- as.tfa_metadata(meta)
#' coverage <- meta$Coverage
#' geno <- coverage_adjust(genotype, coverage, K = 4, log = TRUE)
#'
#' mod  <- tfa(age,
#'             geno,
#'             k = 3,
#'             lambda = 5e-1,
#'             center = TRUE,
#'             coverage = coverage,
#'             log = TRUE)
#'
#' target = c("England_Bell_Beaker", "England_BA")
#' source = c("Anatolia_N", "Russia_Yamnaya", "Serbia_HG")
#' ancestry_coefficients(model = mod,
#'                       metadata = metadata,
#'                       source = source,
#'                       target = target)
#'
#' source = c("Scotland_N", "Russia_Yamnaya", "Serbia_HG")
#' ancestry_coefficients(model = mod,
#'                       metadata = metadata,
#'                       source = source,
#'                       target = target)
#' detach(England_BA)
#' # rm(list = ls())
#' @references François, O., Jay, F. (2020). Factor analysis of ancient DNA samples. Under review.
#' @author Olivier Francois, \email{olivier.francois@@univ-grenoble-alpes.fr}
#' @seealso \code{\link{England_BA}}, \code{\link{tfa}}, \code{\link{coverage_adjust}}, \code{\link{choose_lambda}}
ancestry_coefficients <- function(model, metadata, source, target, individual = FALSE){

  if (class(model) != "tfa") {
    stop("Argument 'model' not of class 'tfa'.")
  }

  if (class(metadata) != "tfa_metadata") {
    stop("Argument 'metadata' not of class 'tfa_metadata'.")
  }

  if (!is.character(source)) {
    stop("Argument 'source' not of class 'character'.")
  }

  if (!is.character(target)) {
    stop("Argument 'target' not of class 'character'.")
  }

  group_ID <- unique(metadata$Group.ID)

  if (sum(target %in% group_ID) != length(target)) {
    stop("a target population was not found. Check target group IDs.")
  }

  if (sum(source %in% group_ID) != length(source)) {
    stop("a source population was not found. Check source group IDs")
  }

  k = length(source)

  if (2 > k){
    stop("At least two sources required.")
  }

  if (k > ncol(model$u) + 1){
    stop("Number of sources greater than the number of factors + 1 in 'model'.")
  } else {
    if (k == ncol(model$u) + 1){
      warning("Number of sources equal to number of factor + 1 may provide inconsistent estimates.")
    }
  }


  ancestry <- NULL

  for (targ in target){

    if (k > 2){
     M <- matrix(NA,  nrow = k - 1, ncol = k)

      for (j in 1:k){
       if (sum(metadata$Group.ID == source[j]) > 1){
         M[,j] <- apply(model$u[metadata$Group.ID == source[j], 1:(k-1)], 2, mean)
       } else {
         M[,j] <- model$u[metadata$Group.ID == source[j], 1:(k-1)]
       }
      }

     if (!individual){

      if (sum(metadata$Group.ID == targ) > 1){
        Z <- apply(model$u[metadata$Group.ID == targ, 1:(k-1)], 2, mean)
      } else {
        Z <- model$u[metadata$Group.ID == targ, 1:(k-1)]
      }

      Z <- t(as.matrix(Z))
      A <- M - M[,1]
      a <- solve(A[,-1]) %*% t(Z - M[,1])
      a <- t(a)
      a <- cbind(1 - sum(a), a)

     } else {
       Z_target <- model$u[metadata$Group.ID == targ, 1:(k-1)]
       nt <- nrow(Z_target)
       a_target <- NULL
       A <- M - M[,1]
       for (i in 1:nt){
         Z <- t(as.matrix(Z_target[i,]))
         ai <- solve(A[,-1]) %*% t(Z - M[,1])
         ai <- t(ai)
         ai <- cbind(1 - sum(ai), ai)
         a_target <- rbind(a_target, ai)
       }
       a <- a_target
     }
    } else { #if nb of sources = 2

      M <- matrix(NA,  nrow = k, ncol = k) # columns are group means (2 groups)
      for (j in 1:k){
        if (sum(metadata$Group.ID == source[j]) > 1){
          M[,j] <- apply(model$u[metadata$Group.ID == source[j], 1:k], 2,
                         mean)
        } else {
          M[,j] <- model$u[metadata$Group.ID == source[j], 1:k]
        }
      }

      if (!individual){
      if (sum(metadata$Group.ID == targ) > 1){
        Z <- apply(model$u[metadata$Group.ID == targ, 1:k], 2, mean)
      } else {
        Z <- model$u[metadata$Group.ID == targ, 1:k]
      }
      Z <- t(as.matrix(Z))

      # computes coordinates in vector basis
      A <- M - M[,1]
      a <- (Z - M[,1]) %*% A[,-1]
      a <- a / t(A[,-1]) %*% A[,-1]
      a <- cbind(1 - a, a)
      } else {
        ## individuals / bug needs correction /done
        Z_target <- model$u[metadata$Group.ID == targ, 1:k] #factors 1:2 for target indiv
        nt <- nrow(Z_target) #number of target indiv
        Z_target <- as.matrix(Z_target)
        a_target <- NULL

        A <- M - M[,1]

        for (i in 1:nt){
          ai <- (Z_target[i,] - M[,1]) %*% A[,-1]
          ai <- ai / t(A[,-1]) %*% A[,-1]
          ai <- cbind(1 - sum(ai), ai)
          a_target <- rbind(a_target, ai)
        }
        a <- a_target
       }
    }
    ancestry <- rbind(ancestry, a)
  }

  colnames(ancestry) <- source

  if (!individual){
    rownames(ancestry) <- target
  } else {
    rownames(ancestry) <- metadata$Instance.ID[metadata$Group.ID %in% target]
  }

  return(ancestry)
}
