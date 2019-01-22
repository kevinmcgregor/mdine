#' Microbiome family-level data from Crohn's and control samples
#'
#' Counts from 16S sequencing data on treatment-naive subjects with Crohn's disease
#' or controls.  Other covariates available include age, sex, and antibiotic use.
#' 
#' @docType data
#'
#' @usage data(crohns)
#'
#' @format Contains a dataframe called "covars" which contains the covariate data, as well as a matrix called "otu.counts" which contains the family-level counts of each subject.
#'
#' @keywords datasets crohn's microbiome
#'
#' @references Gevers et al. (2014) Cell Host & Microbe 15.3 (2014): 382-392.
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/24629344}{PubMed})
#'
#' @source \href{https://qiita.ucsd.edu/study/description/1939}{Qiita database}
#'
#' @examples
#' data(grav)
#' 
"crohns"