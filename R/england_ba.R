#' Ancient human DNA samples from David Reich's lab
#' @description {The data include 10,000 filtered SNP genotypes from 118 samples, from Neolithic and Bronze Age Great Britain,
#' Steppe (Yamnaya), early farmers from Anatolia, and hunter-gatherers from Serbia. The data set is structured as a list
#' named 'England_BA'
#'
#'$age: vector of ages in  years cal BP.
#'
#'$genotype: a numeric matrix containing 118 genotypes. The SNPs were randomly selected from a much larger set of variants.
#'
#'$meta: metadata on individual ID, group ID, country of origin and coverage.
#'
#'}
#' @name england_ba
#' @docType data
#'
NULL



#' Ancient human DNA samples from David Reich's lab
#'
#' A dataset containing 10,000 filtered SNP genotypes from 118 samples, from Neolithic and Bronze Age Great Britain,
#' Steppe (Yamnaya), early farmers from Anatolia, and hunter-gatherers from Serbia. Metadata and sample
#' ages are also provided.
#'
#' @format List with 3 objects
#' \describe{
#'   \item{age}{ages, in years cal BP}
#'   \item{genotype}{a numeric matrix containing diploid genotypes (0,1,2) for
#'   118 ancient Eurasians. The SNPs were randomly selected from a much larger
#'   set of polymorphisms}
#'   \item{meta}{metadata on individuals: ID, group ID, country of origin and coverage}
#' }
#' @source \url{https://reich.hms.harvard.edu/datasets}
"England_BA"
