#' Factor analysis of temporal or ancient DNA samples
#'
#' \code{tfa} estimates factors describing population structure for temporal samples
#' of DNA, correcting individual scores for the effect of allele frequency drift through time
#'
#' @param sample_ages a numeric vector corresponding to the age of each sample where age = 0 is for
#' present-day individuals. By default, ages are converted into normalized dates between 0 and 1
#' (date = 1 for present-day individuals).
#' @param Y an nxp numeric matrix containing genetic information for n individuals recorded in p columns.
#' Genetic information could be encoded by any numeric value, not necessarily an integer value. Missing data are not allowed.
#' @param lambda a nonnegative numeric value which corresponds to the drift parameter (noise-to-temporal-signal ratio).
#' @param k an integer value for the number of factor to compute. The default value is k = 2.
#' @param cov_matrix a user-specified prior covariance matrix for the n samples. If \code{NULL},
#' the prior matrix is computed as
#' \code{C[i,j] = min(t_i, t_j)},
#' where \code{t_i} and \code{t_j} are the normalized sample dates, taking values between 0 and 1. The option is
#' useful when the sample size is large, and when the user wants to use a pre-computed covariance matrix
#' to save time.
#' @param center a logical value indicating whether the genetic values should be centered by substracting
#' their row mean.
#' @param coverage a numerical vector containing information on DNA sample coverage. Note that
#' coverage differences might strongly bias factor analysis results. Including coverage information
#' allows the program to adjust for coverage bias by using a local regression (\code{loess}).
#' We also suggest that correction of the data for low coverage should be performed before analysis by using the function \code{coverage_adjust}.
#' @param log a logical value indicating that corrections are performed from \code{log(coverage)}
#'  instead of \code{coverage}.
#'
#' @return A list with the following attributes:
#' \describe{
#'    \item{u}{an nxk numeric matrix containing the k corrected factors}
#'    \item{singular.values}{a vector of size n containing the singular values (std dev) for the adjusted factors}
#'    \item{tn}{a vector containing the normalized sample dates}
#'    \item{cov}{the covariance matrix C used for correction (default: Brownian covariance matrix)}
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
#' legend(x = 10, y = -19, cex = .6,
#'       legend = c("Early Farmers", "Hunter Gatherers", "Steppe"),
#'       col = c("salmon3", "olivedrab", "darkblue"), pch = 8)
#' legend(x = -24, y = -19, cex = .6,
#'        legend = c("Neolithic GB", "Bronze Age GB", "Bell Beaker"),
#'        col = c("salmon1", "yellow3", "yellow4"), pch = 19)
#' detach(England_BA)
#' # rm(list = ls())
#' @references François, O., Liégeois, S., Demaille, B., Jay, F. (2019). Inference of population genetic structure from temporal samples
#' of DNA. bioRxiv, 801324. \url{https://www.biorxiv.org/content/10.1101/801324v3}
#'
#' François, O., Jay, F. (2020). Factor analysis of ancient DNA samples. Under review.
#'
#' @seealso \code{\link{England_BA}}, \code{\link{coverage_adjust}}, \code{\link{choose_lambda}}, \code{\link{ancestry_coefficients}}
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
  if (nrow(Y) != n) stop("Number of rows (samples) of Y not equal to the length of
                           sample ages.")


  # Covariance model
  # Default is Brownian covariance matrix
  if (is.null(cov_matrix)){

  # conversion of ages as times: tn in (0,1)
    range_ages <- max(sample_ages) - min(sample_ages)
    tn <- 1 - (sample_ages - min(sample_ages))/range_ages

  # sample time normalisation
    var_Y <- apply(Y, 1, FUN = var)
    tn <- min(var_Y) +  (max(var_Y) - min(var_Y)) *tn

  # Brownian covariance model for scaled data
    C <- matrix(NA, n, n)
      for (i in 1:n){
        for (j in 1:n)  C[i,j] <- min(tn[i], tn[j])}
      } else {
  # Use cov_matrix as covariance matrix
    C <- as.matrix(cov_matrix)
     if (nrow(C) != n | ncol(C) != n){
      stop("'cov_matrix' has incorrect dimensions.")
      }
 # conversion of ages
    range_ages <- max(sample_ages) - min(sample_ages)
    tn <- 1 - (sample_ages - min(sample_ages))/range_ages
 # normalisation
    var_Y <- apply(Y, 1, FUN = var)
    tn <- min(var_Y) +  (max(var_Y) - min(var_Y)) *tn
      }

 # scale Y
  if (center) {
   Y = t(scale(t(Y), scale = FALSE, center = TRUE))
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


#' Graphical heuristic for choosing the drift parameter in \code{tfa()}
#'
#' \code{choose_lambda} runs temporal factor analyses for drift parameters in a specified
#' range, and ouputs the percentage of variance of sample time explained by all factors for each value.
#'
#' @param model an object of class 'tfa' adjusted with \code{tfa()}.
#' @param Y an nxp numeric matrix containing genetic information for n individuals recorded in p columns.
#' Genetic information could be encoded by any numeric value, not necessarily an integer value.
#' Missing data are not allowed.
#' @param min_range a numeric value for the minimal range of drift parameter tested. Log 10 scale.
#' The default value is \code{min_range = -5}.
#' @param max_range a numeric value for the maximal range of drift parameter tested. Log 10 scale.
#' The default value is \code{max_range = 0.1}.
#' @param grid_size an integer value for the number of drift parameter tested. The default is
#' \code{grid_size = 10}.
#' @param plot_res a logical indicating whether the results should be displayed graphically or not. The
#' default value is \code{TRUE}.
#' @return A vector of percentages of variance of sample time explained by factors for each value of
#' the drift parameter in the specified range.
#' @details This function requires that a preliminary model has been adjusted with K factors, where K is
#' the number of source populations. The curve exhibits a drop for some value in the range, and the heuristic
#' consists of choosing the largest values corresponding to removing time variation in factor K.
#' @author Olivier Francois, \email{olivier.francois@@univ-grenoble-alpes.fr}
#' @importFrom stats lm
#' @export
#' @examples
#' library(tfa)
#'
#' # Ancient DNA from Bronze Age Great Britain samples
#' data(england_ba)
#' attach(England_BA)
#'
#' coverage <- meta$Coverage
#' geno <- coverage_adjust(genotype, coverage, K = 4, log = TRUE)
#'
#' # Adjust an FA model
#' mod  <- tfa(age,
#'             geno,
#'             k = 3,
#'             lambda = 5e-1,
#'             center = TRUE,
#'             coverage = coverage,
#'             log = TRUE)
#'
#' r_2 <- choose_lambda(mod,
#'                      geno,
#'                      min_range = -4,
#'                      max_range = 3)
#' abline(v=log10(5e-1), col = "orange", lty = 2)
#'
#' # Remove HG from Serbia to keep k = 2 ancestral populations
#'
#' age2 <- age[meta$Group.ID != "Serbia_HG"]
#' geno2 <- geno[meta$Group.ID != "Serbia_HG",]
#'
#' # Adjust an FA model
#' mod  <- tfa(age2,
#'             geno2,
#'             k = 2,
#'             lambda = 5e-2,
#'             center = TRUE)
#'
#' r_2 <- choose_lambda(mod,
#'                      geno2,
#'                      min_range = -3,
#'                      max_range = 2)
#' abline(v=log10(5e-2), col = "orange", lty = 2)
#' detach(England_BA)
#' # rm(list = ls())
#' @references François, O., Liégeois, S., Demaille, B., Jay, F. (2019). Inference of population genetic structure from temporal samples
#' of DNA. bioRxiv, 801324. \url{https://www.biorxiv.org/content/10.1101/801324v3}
#' @seealso \code{\link{England_BA}}, \code{\link{tfa}}
choose_lambda <- function(model,
                          Y,
                          min_range = -5,
                          max_range =  1,
                          grid_size = 10,
                          plot_res = TRUE){

  if(!inherits(model, what = "tfa")){
    stop("Object 'model' not of class 'tfa'.")
  }

  tn <- model$adjusted.times
  n <- length(tn)
  if (!is.matrix(Y)) stop("Y must be a matrix.")
  if (anyNA(Y)) stop("NA not allowed in genetic data.")
  if (nrow(Y) != n) stop("Number of rows (samples) of Y not equal to the length of
                         sample ages.")



  r_squared <- NULL

  # define a grid of lambda values on which FA models are adjusted
   llambda <- seq(min_range, max_range, length = grid_size)
   lambda <- 10^llambda

  # times and covariance
   C <- model$cov
   k <- ncol(model$u)


  for (la in lambda){
        fa <- tfa(sample_ages = 1 - tn,
                     Y,
                     k = k,
                     lambda = la,
                     cov_matrix = C,
                     center = TRUE)

    mod_lm <- lm(tn ~ fa$u)
    r_squared <- c(r_squared, summary(mod_lm)$adj.r.squared)
  }
  if (plot_res){
    plot(llambda,
           r_squared,
           xlab = "Log parameter",
           ylab = "Var. explained",
           cex = 1)
    points(llambda,
           r_squared,
           col = "lightblue",
           pch = 19,
           cex = .8)
  }
  return(r_squared)
}




#' Compute proportions of variance explained by temporal
#' distortion at each locus
#'
#' \code{prop_variance} returns a vector of squared correlation coefficients representing the
#' proportion of variation captured by the specified covariance model at each locus.
#' @param model an object of class 'tfa' adjusted with \code{tfa()}.
#' @param Y an nxp numeric matrix containing genetic information for n individuals recorded in p columns.
#' Genetic information could be encoded by any numeric value, not necessarily an integer value.
#' Missing data are not allowed. The matrix should have the same number of rows as in the adjusted model.
#' @return A vector of probabilities representing the proportion of variation captured the specified
#' covariance model at each locus.
#' @details This function requires that a preliminary model has been adjusted with K factors, where K is
#' the number of source populations.
#' @author Olivier Francois, \email{olivier.francois@@univ-grenoble-alpes.fr}
#' @importFrom stats lm
#' @importFrom graphics plot points
#' @export
#' @examples
#' library(tfa)
#'
#' # Ancient DNA from Bronze Age Great Britain samples
#' data(england_ba)
#'
#' attach(England_BA)
#'
#' # Remove HG from Serbia to keep k = 2 ancestral populations
#' age <- age[meta$Group.ID != "Serbia_HG"]
#' geno <- genotype[meta$Group.ID != "Serbia_HG",]
#'
#' # Adjust an FA model
#' mod  <- tfa(age,
#'             geno,
#'             k = 2,
#'             lambda = 1e-3,
#'             center = TRUE)
#'
#' r_squared <- prop_variance(mod, geno)
#' summary(r_squared)
#' hist(r_squared, col = "darkviolet")
#' detach(England_BA)
#' # rm(list = ls())
#' @references François, O., Liégeois, S., Demaille, B., Jay, F. (2019). Inference of population genetic structure from temporal samples
#' of DNA. bioRxiv, 801324. \url{https://www.biorxiv.org/content/10.1101/801324v3}
#' @seealso \code{\link{England_BA}}, \code{\link{tfa}}
prop_variance <- function(model,
                          Y){

  if(!inherits(model, what = "tfa")){
    stop("Object 'model' not of class 'tfa'.")
  }

  n <- length(model$adjusted.times)
  if (!is.matrix(Y)) stop("Y must be a matrix.")
  if (anyNA(Y)) stop("NA not allowed in genetic data.")
  if (nrow(Y) != n) stop("Number of rows (samples) of Y not equal to the length of
                         sample ages.")

  ## get covariance model from the adjusted object
  C <- model$cov
  eg <- eigen(C)

  ## reduce the number of eigenvectors
  k1 <- min(which(cumsum(eg$values)/sum(eg$values) > 0.995))
  k <- ncol(model$u)

  v <- apply(Y, 2, var)
  if (sum(v==0)>0){
    warning("Columns of the Y matrix have null variance.")
  }

  ## variance explained by the FA model
  mod_lm <- lm(Y ~ ., data = data.frame(eg$vectors[,1:k1], model$u[,1:k]))
  sm <- summary(mod_lm)
  v_all <- sapply(sm, FUN = function(x) x$r.squared)

  ## variance explained by the k factors
  mod_lm <- lm(Y ~ ., data = data.frame(model$u[,1:k]))
  sm <- summary(mod_lm)
  v_fa <- sapply(sm, FUN = function(x) x$r.squared)

  ## residual variance (explained by the eigenvectors)
  r_squared <- v_all - v_fa
  return(r_squared)
}

