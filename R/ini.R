#' @importFrom stats as.formula cor rbeta rnorm terms
NULL

#' Sample Genotype C17
#'
#' Chromosome 17 region  q12 taken from 1000 genome project,  with 10146 common
#' SNP (MAF  > 0.05), no  missings, and 1870  unrelated individuals out  of the
#' original 2502.
#'
#' The sample data is organized into a 1870 x 10146 matrix.
#' 
#' @name C17
#' @docType data
#' @author Xiaoran Tong \email{xiaoran.tong.cn@gmail.com}
#' @references
#' \url{ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502}
#' @keywords data
NULL

#' Sample Allele Frequency
#'
#' The minor allele frequency of sample genotype data -- "C17".
#' 
#' @name F17
#' @docType data
#' @author Xiaoran Tong \email{xiaoran.tong.cn@gmail.com}
#' @keywords data
#' @seealso [C17]
NULL

utils::globalVariables(c("C17", "F17"))
