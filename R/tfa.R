#' Factor analysis of temporal and ancient individual genotypes
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
#'    \item{tn}{vector of normalized sample dates}
#'    \item{u}{an nxk numeric matrix containing the k corrected factors}
#'    \item{cov}{the Brownian covariance matrix used for correction}
#' }
#'
#' @export
#' @examples
#' library(tfa)
#'
#' # Ancient DNA from Bronze Age Great Britain samples
#' # including Steppe (Yamnaya) and early farmers from Anatolia
#' data(england_ba)
#' attach(England_BA)
#' coverage <- as.numeric(as.character(meta$Coverage))
#' mod  <- tfa(age,
#'            genotype,
#'            k = 2,
#'            lambda = 1e-3,
#'            center = TRUE,
#'            coverage = coverage,
#'            log = TRUE)
#'
#' plot(mod$u, pch = 19, cex = 2, col = "grey90",
#'      xlab = "Factor 1", ylab = "Factor 2",
#'      ylim = c(-26,20))
#'
#' points(mod$u[meta$Country == "Great Britain",],
#'        pch = 19, cex = 1, col = "yellow3")
#' points(mod$u[meta$Country == "Russia",],
#'       pch = 8, cex = .6, col = "darkblue")
#' points(mod$u[meta$Country == "Turkey",],
#'       pch = 8, cex = .6, col = "salmon3")
#' points(mod$u[meta$Group.ID == "England_MBA",],
#'       pch = 8, cex = .6, col = "yellow4")
#' abline(h = 0, lty = 2, col = "orange")
#' legend(x = -55, y = -17,
#'      c("Yamnaya","England_EBA  ","England_MBA   ","Anatolia"),
#'      horiz = TRUE,
#'      col = c("darkblue", "yellow3", "yellow4", "salmon3"),
#'      pch = c(8,19,8,8), cex = .8, bty = "n")
#'
#' ne = mean(mod$u[meta$Country == "Turkey",1])
#' ru = mean(mod$u[meta$Country == "Russia",1])
#'
#' # Fraction of Yamnaya ancestry in England Bronze Age genomes
#' gb = mean(mod$u[meta$Country == "Great Britain",1])
#' (ne - gb)/(ne - ru)
#'
#' # Fraction of Yamnaya ancestry in Middle Bronze Age genomes
#' mba = mean(mod$u[meta$Group.ID == "England_MBA",1])
#' (ne - mba)/(ne - ru)
#'
#' detach(England_BA)
#' @references François, O., Liégeois, S., Demaille, B., Jay, F. (2019). Inference of population genetic structure from temporal samples
#' of DNA. bioRxiv, 801324. \url{https://www.biorxiv.org/content/10.1101/801324v3}
#' @seealso \code{\link{england_ba}}
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
        #plot(log(coverage), Wn[,k], col = "grey", pch = 19, cex = .8)
        mod = loess(Un[,k]  ~ log(coverage))
        #points(log(coverage), fitted(mod), col = 2, pch = 19, cex = 1)
        Un[,k] <- Un[,k] - fitted(mod)}
      Un <-  prcomp(Un)$x
    } else {
      for (i in 1:k){
        #plot(coverage, Wn[,k], col = "grey", pch = 19, cex = .8)
        mod = loess(Un[,k]  ~ coverage)
        #points(coverage, fitted(mod), col = 2, pch = 19, cex = 1)
        Un[,k] <- Un[,k] - fitted(mod)}
      Un <-  prcomp(Un)$x
    }
  }
  # returns corrected factors in $u
  return(list(tn = tn, u = Un, cov = C))
}
