
#' Get adjacency matrix for estimated networks based on credible intervals
#'
#' @param obj An object of 'mdine' class
#' @param weighted Logical: should the adjacency matrices be weighted according to the corresponding partial correlation matrix?
#'
#' @return A list containing adjacency matrices for the two estimated precision matrices
#' @export
#'
#' @examples ls()
ci2adj <- function(obj, weighted=FALSE) {
  if (class(obj) != "mdine") stop("obj must be of class \"mdine\"")

  adj0 <- (obj$ci$invsigma0[[1]]>0) | (obj$ci$invsigma0[[2]]<0)
  adj1 <- (obj$ci$invsigma1[[1]]>0) | (obj$ci$invsigma1[[2]]<0)
  diag(adj0) <- 0
  diag(adj1) <- 0

  return(list(adj0=adj0, adj1=adj1))
}

