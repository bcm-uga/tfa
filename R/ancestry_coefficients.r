#' Convert a data frame as an object of class tfa_data
#'
#' \code{as.tfa_data} checks the contents of a metadata file for the presence of group IDs, and
#' convert the data frame as an object of class tfa_data.
#'
#' @param dataframe a dataframe containing metadata for a factor analysis of ancient DNA. The meta
#' data could contain information on sample IDs, Group.IDs, country of origin, sample age, etc. Only
#' the column with name 'Group.ID' is mandatory.
#' @return an object of class tfa_data
#' @export
#' @examples
#' library(tfa)
#'
#' # Ancient DNA from Bronze Age Great Britain samples
#' # including Steppe (Yamnaya), hunter gatherers, and early farmers from Anatolia
#' data(england_ba)
#' attach(England_BA)
#' metadata <- as.tfa_data(meta)
#' detach(England_BA)
#' @references François, O., Jay, F. (2020). Factor analysis of ancient DNA samples.
#' @seealso \code{\link{England_BA}}, \code{\link{ancestry_coefficients}}
as.tfa_data <- function(dataframe){

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
  obj <- dataframe
  class(obj) <- "tfa_data"
  return(obj)
}

#' Compute ancestry coefficients from sources populations for a target group
#'
#' \code{ancestry_coefficients} compute ancestry coefficients from sources populations for a target group
#'
#' @param model an object of class tfa with the same number of individuals as in metadata.
#' @param metadata an object of class tfa_data containing metadata on Group.IDs.
#' @param source a vector of (two or three) ancestral group IDs (character).
#' @param target a target group ID.
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
#' metadata <- as.tfa_data(meta)
#' coverage <- meta$Coverage
#' geno <- coverage_adjust(genotype, coverage, K = 3, log = TRUE)
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
#' @references François, O., Jay, F. (2020). Factor analysis of ancient DNA samples.
#' @author Olivier Francois, \email{olivier.francois@@univ-grenoble-alpes.fr}
#' @seealso \code{\link{England_BA}}, \code{\link{tfa}}
ancestry_coefficients <- function(model, metadata, source, target){

  if (class(model) != "tfa") {
    stop("Argument 'model' not of class 'tfa'.")
  }

  if (class(metadata) != "tfa_data") {
    stop("Argument 'metadata' not of class 'tfa_data'.")
  }

  if (!is.character(source)) {
    stop("Argument 'source' not of class 'character'.")
  }

  if (!is.character(target)) {
    stop("Argument 'target' not of class 'character'.")
  }

  group_ID <- unique(metadata$Group.ID)

  if (sum(target %in% group_ID) != length(target)) {
    stop("'target' population not found.")
  }

  k = length(source)

  if (2 > k){
    stop("At least two sources required.")
  }

  if (k > ncol(model$u) + 1){
    stop("Number of sources greater than the number of factors + 1 in 'model'.")
  }

  ancestry <- NULL

  for (targ in target){

    if (k > 2){
     M <- matrix(NA,  nrow = k - 1, ncol = k)

      for (j in 1:k){
       M[,j] <- apply(model$u[metadata$Group.ID == source[j], 1:(k-1)], 2, mean)
      }

      Z <- apply(model$u[metadata$Group.ID == targ, 1:(k-1)], 2, mean)
      Z <- t(as.matrix(Z))

      A <- M - M[,1]
      a <- solve(A[,-1]) %*% t(Z - M[,1])
      a <- t(a)
      a <- cbind(1 - sum(a), a)
    } else {

      M <- matrix(NA,  nrow = k, ncol = k)
      for (j in 1:k){
       M[,j] <- apply(model$u[metadata$Group.ID == source[j], 1:k], 2, mean)
     }

      Z <- apply(model$u[metadata$Group.ID == targ, 1:k], 2, mean)
      Z <- t(as.matrix(Z))
      A <- M - M[,1]
      a <- (Z - M[,1]) %*% A[,-1]
      a <- a / t(A[,-1]) %*% A[,-1]
      a <- cbind(1 - a, a)
    }
    ancestry <- rbind(ancestry, a)
  }

  colnames(ancestry) <- source
  rownames(ancestry) <- target
  return(ancestry)
}