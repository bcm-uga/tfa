#' Factor analysis of temporal or ancient DNA samples
#'
#' \code{tfa} estimates factors describing population structure for temporal samples
#' of DNA, correcting individual scores for the effect of allele frequency drift through time
#'
#' @param sample_ages a numeric vector corresponding to the ages of each sample where age = 0 for
#' present-day individuals). By default, ages are converted into dates between 0 and 1
#' (date = 1 for present-day individuals).
#' @param Y an nxp numeric matrix containing genetic information for n individuals recorded in p columns.
#' Genetic information could be encoded as any numeric value, not necessarily an integer value. Missing data are not allowed.
#' @param lambda a nonnegative numeric value which corresponds to a scale parameter (noise-to-temporal-signal ratio).
#' @param k an integer value for the number of factor to compute. The default value is k = 2.
#' @param cov_matrix a user-specified nxn prior covariance matrix for the n samples. If \code{NULL},
#' the prior matrix is computed as
#' \code{C[i,j] = min(t_i, t_j)},
#' where t_i and t_j are the normalized sampled dates, taking values between 0 and 1. The option is
#' useful when the sample size is large, and when the user wants to avoid recomputing a same C
#' matrix in every run.
#' @param center a logical value indicating whether the genetic values should be centered by substracting
#' their row mean.
#' @param coverage a numerical vector containing information on DNA sample coverage. Note that
#' coverage differences might strongly bias factor analysis results. Including coverage information
#' allows the program to correct coverage bias on factors via local regression (\code{loess}).
#' We also suggest that correction of the data for low coverage should be performed before analysis
#' with \code{tfa}, for example, using the function \code{coverage_adjust}.
#' @param log a logical value indicating that corrections are performed from log(coverage) instead of coverage.
#'
#' @return A list with the following attributes:
#' \describe{
#'    \item{u}{an nxk numeric matrix containing the k corrected factors}
#'    \item{singular.values}{a vector of size n containing the singular values}
#'    \item{tn}{vector of normalized sample dates}
#'    \item{cov}{the Brownian covariance matrix used for correction}
#' }
#' @author Olivier Francois, \email{olivier.francois@@univ-grenoble-alpes.fr}
#' @importFrom stats prcomp loess fitted var
#' @export
#' @examples
#' library(tfa)
#'
#' # Ancient DNA from Bronze Age Great Britain samples
#' # including Yamnaya, early farmers (Anatolia) and hunter-gatherers (Serbia)
#' data(england_ba)
#'
#' attach(England_BA)
#' coverage <- meta$Coverage
#' geno <- coverage_adjust(genotype, coverage, K = 3, log = TRUE)
#'
#' mod  <- tfa(age,
#'             geno,
#'             k = 2,
#'             lambda = 5e-1,
#'             center = TRUE,
#'             coverage = coverage,
#'             log = TRUE)
#'
#' plot(mod$u, pch = 19, cex = 2, col = "grey90",
#'      xlab = "Factor 1", ylab = "Factor 2",
#'      main = "FA")
#'
#' m_yamnaya <- apply(mod$u[meta$Group.ID == "Russia_Yamnaya",],
#'                    2, mean)
#' m_anatolia <- apply(mod$u[meta$Group.ID == "Anatolia_N",],
#'                     2, mean)
#' m_hg <- apply(mod$u[meta$Group.ID == "Serbia_HG",],
#'               2, mean)
#'
#' points(rbind(m_yamnaya, m_anatolia, m_hg), lwd = 2)
#'
#' lines(rbind(m_yamnaya, m_anatolia, m_hg, m_yamnaya))
#'
#' points(mod$u[meta$Group.ID == "Russia_Yamnaya",],
#'        pch = 8, cex = .6, col = "darkblue")
#'
#' points(mod$u[meta$Group.ID == "Anatolia_N",],
#'        pch = 8, cex = .6, col = "salmon3")
#'
#' points(mod$u[meta$Group.ID ==  "Serbia_HG",],
#'        pch = 8, cex = .6, col = "olivedrab")
#'
#' points(mod$u[meta$Group.ID == "England_Bell_Beaker",],
#'        pch = 19, cex = .6, col = "yellow4")
#'
#' points(mod$u[meta$Group.ID == "England_BA",],
#'        pch = 19, cex = .6, col = "yellow3")
#'
#' points(mod$u[meta$Group.ID %in% c("England_N", "Scotland_N"),],
#'        pch = 19, cex = .6, col = "salmon1")
#'
#' legend(x = 10, y = -15, cex = .6,
#'        legend = c("Early Farmers", "Hunter Gatherers", "Steppe"),
#'        col = c("salmon3", "olivedrab", "darkblue"), pch = 8)
#'
#' legend(x = 10, y = 20.5, cex = .6,
#'        legend = c("Neolithic GB", "Bronze Age GB", "Bell Beaker"),
#'        col = c("salmon1", "yellow3", "yellow4"), pch = 19)
#' detach(England_BA)
#' # rm(list = ls())
#' @references François, O., Liégeois, S., Demaille, B., Jay, F. (2019). Inference of population genetic structure from temporal samples
#' of DNA. bioRxiv, 801324. \url{https://www.biorxiv.org/content/10.1101/801324v3}
#' @seealso \code{\link{England_BA}}, \code{\link{coverage_adjust}}
tfa <- function(sample_ages,
                Y,
                k = 2,
                lambda = 1e-3,
                cov_matrix = NULL,
                center = TRUE,
                coverage = NULL,
                log = TRUE){

  n <- length(sample_ages)
  if (!is.matrix(Y)) stop("Y must be a matrix.")
  if (anyNA(Y)) stop("NA not allowed in genetic data.")
  if (anyNA(sample_ages)) stop("NA not allowed in sample ages.")
  if (dim(Y)[1] != n) stop("Number of Y rows (samples) not equal to the length of
                           sample ages.")


  # Covariance model
  # Default is Brownian covariance matrix
  if (is.null(cov_matrix)){

  # conversion of ages as tn in (0,1)
  range_ages <- max(sample_ages) - min(sample_ages)
  tn <- 1 - (sample_ages - min(sample_ages))/range_ages
  # normalisation of times
  var_Y <- apply(Y, 1, FUN = var)
  tn <- min(var_Y) +  (max(var_Y) - min(var_Y)) *tn

  # scale Y
  if (center) {
    Y = t(scale(t(Y), scale = FALSE, center = TRUE))
  }
  # Brownian model for scaled data
  C <- matrix(NA, n, n)
  for (i in 1:n){
    for (j in 1:n)  C[i,j] <- min(tn[i], tn[j])}
  } else
  {
    C <- as.matrix(cov_matrix)
    if (dim(C) != c(n,n)){
      stop("Covariance matrix in 'cov_matrix' has incorrect dimensions.")
    }
  }

  ### covariance eigenspectrum
  ec <- eigen(C)
  Qn <- ec$vect
  lambda_n <- ec$values

  # correction computation
  D <- diag( sqrt(lambda/(lambda_n + lambda)) )
  D_inv <- diag( sqrt((lambda_n + lambda)/lambda) )

  sv <- svd(D %*% t(Qn) %*% Y, nu = k)
  Un <- Qn %*% D_inv %*% sv$u %*% diag(sv$d[1:k])
  sw <- svd(Un %*% t(sv$v[,1:k]), nu = k)
  Un <- sw$u %*% diag(sw$d[1:k])

  # Correction for coverage using loess regression
  if (!is.null(coverage)){
    if (log){
      for (i in 1:k){
        mod = loess(Un[,k]  ~ log(coverage))
        Un[,k] <- Un[,k] - fitted(mod)}
      Un <-  prcomp(Un)$x
    } else {
      for (i in 1:k){
        mod = loess(Un[,k]  ~ coverage)
        Un[,k] <- Un[,k] - fitted(mod)}
      Un <-  prcomp(Un)$x
    }
  }
  # returns corrected factors in $u
  obj <- list(u = Un, singular.values = sw$d, adjusted.times = tn, cov = C)
  class(obj) <- "tfa"
  return(obj)
}
