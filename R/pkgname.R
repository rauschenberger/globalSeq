
#' @title
#' globalSeq
#' 
#' @description
#' Testing for association between
#' RNA-Seq and high-dimensional data
#' 
#' \strong{Background}
#' \cr Testing for association between RNA-Seq
#' and other genomic data is challenging
#' due to high variability of the former
#' and high dimensionality of the latter.
#' 
#' \strong{Results}
#' \cr Using the negative binomial distribution and a random effects model,
#' we developed an omnibus test that overcomes both difficulties.
#' It may be conceptualised as a test of overall significance
#' in regression analysis, where the response variable is overdispersed
#' and the number of explanatory variables exceeds the sample size.
#' 
#' \strong{Conclusions}
#' \cr The proposed method can detect genetic and epigenetic alterations
#' that affect gene expression. It can examine complex regulatory
#' mechanisms of gene expression.
#' 
#' \strong{Reference}
#' \cr A Rauschenberger, MA Jonker, MA van de Wiel, and RX Menezes (2016).
#' "Testing for association between RNA-Seq and high-dimensional data",
#' \emph{BMC Bioinformatics}. 17:118.
#' \href{http://dx.doi.org/10.1186/s12859-016-0961-5}{html}
#' \href{http://www.biomedcentral.com/content/pdf/s12859-016-0961-5.pdf}{pdf}
#' (open access)
#' 
#' @seealso
#' The following command opens the vignette:
#' \cr \code{utils::vignette("globalSeq")}
#' \cr
#' \cr \code{\link{omnibus}} tests entire covariate sets
#' \cr \code{\link{proprius}} shows individual contributions
#' \cr \code{\link{cursus}} analyses the whole genome
#'
#' @keywords documentation
#' @docType package
#' @name globalSeq
#' @aliases globalseq
NULL 
