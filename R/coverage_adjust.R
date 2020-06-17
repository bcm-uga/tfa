#' Adjust ancient genotypes for coverage bias
#'
#' \code{coverage_adjust} adjusts the matrix of genotypes for low and uneven coverage of
#' individual genomes.
#'
#' @param Y an nxp numeric data matrix containing genetic information for n individuals recorded in p columns.
#' Genetic information could be encoded as any numeric value, usually as genotypes taking values in 0,1,2.
#' Missing data are not allowed.
#' @param coverage a numerical vector containing information on each DNA sample coverage.
#' @param K a nonnegative integer for the number of components in the data matrix Y.
#' It could be obtained from a PCA of the data matrix.
#' @param lambda  a ridge regularization parameter, to be kept at a small values (defaut 1e-5)).
#' @param log a logical value indicating that corrections are performed from log(coverage) instead of coverage.
#' @return an nxp numeric matrix containing the n corrected genotypes
#' @export
#' @examples
#' library(tfa)
#'
#' # Ancient DNA from Bronze Age Great Britain samples
#' # including Steppe (Yamnaya), hunter gatherers, and early farmers from Anatolia
#' data(england_ba)
#' attach(England_BA)
#' coverage <- meta$Coverage
#' geno <- coverage_adjust(genotype, coverage, K = 4, log = TRUE)
#' detach(England_BA)
#' @references François, O., Jay, F. (2020). Factor analysis of ancient population genomic samples. Under review.
#' @author Olivier Francois, \email{olivier.francois@@univ-grenoble-alpes.fr}
#'
#' @seealso \code{\link{England_BA}}, \code{\link{tfa}}, \code{\link{choose_lambda}}
coverage_adjust <- function(
                  Y = NULL,
                  coverage = NULL,
                  K,
                  lambda = 1e-5,
                  log = FALSE){
    ## Check response matrix Y
    if (is.null(Y)){
      stop("NULL values for the genotype matrix Y.")
    }
    if (!is.matrix(Y)){
      Y <- as.matrix(Y)
    }
    if (anyNA(Y)) {
      stop("The genotype matrix Y contains NA.")
    }

  ## Check coverage values
  X <- coverage
    if (is.null(X)){
      stop("No values for argument coverage.")
    }
    if (!is.numeric(X)){
      stop("Numeric values expected for coverage.")
    }
    if (anyNA(X)) {
      stop("The coverage vector contains NA.")
    }
   if (length(X) != nrow(Y)) {
    stop("Length of coverage is not equal to nrow(Y).")
  }

  n = nrow(Y)

  if (log){
    X <- log(X)
  }

  # run SVD of X and implement Caye's transform
  svx <- svd(x = scale(X, scale = FALSE), nu = n)
  Q <- svx$u

  d_lambda <- c(sqrt(lambda/(lambda + svx$d)), rep(1, n-1))
  D  <- diag(d_lambda)
  d_lambda_inv <- c(sqrt((lambda + svx$d)/lambda), rep(1, n-1))
  D_inv <- diag(d_lambda_inv)

  # run rank K approximation of transformed Y
  svk <- svd(D %*% t(Q) %*% scale(Y, scale = FALSE), nu = K)

  # compute the latent matrix W
  W <- Q %*% D_inv %*% svk$u %*% diag(svk$d[1:K]) %*% t(svk$v[,1:K])

  # output W
  return(W)
  }
