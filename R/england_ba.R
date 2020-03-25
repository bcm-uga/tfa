#' Ancient human DNA samples from David Reich's lab
#'
#' A dataset containing 10,000 filtered SNP genotypes from 118 samples, from Neolithic and Bronze Age Great Britain,
#' Steppe (Yamnaya), early farmers from Anatolia, and hunter-gatherers from Serbia. Metadata and sample
#' ages are also provided.
#'
#' @format List with 3 objects
#' \describe{
#'   \item{age}{a numeric vector containing the sample ages, where age is measured in year cal BP}
#'   \item{genotype}{a numeric matrix containing diploid genotypes (0,1,2) for 118 ancient Eurasians. The 10,000 filtered SNPs were randomly selected from a large set of variants. The matrix contains no missing genotypes.}
#'   \item{meta}{metadata on individuals: \code{Instance.ID}, \code{Group.ID}, \code{Country} of origin and \code{coverage}}
#' }
#' @source \url{https://reich.hms.harvard.edu/datasets}
"England_BA"
