
#' Get adjacency matrix for estimated networks based on credible intervals
#'
#' @param obj An object of class \code{mdine}
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

  if (weighted) {
    adj0 <- -cov2cor(obj$post_mean$invsigma0)*adj0
    adj1 <- -cov2cor(obj$post_mean$invsigma1)*adj1
  }

  diag(adj0) <- 0
  diag(adj1) <- 0

  return(list(adj0=adj0, adj1=adj1))
}


#' Get which elements of the precision matrices differ significantly
#'
#' @param obj An object of class \code{mdine}
#'
#' @return A matrix with a 1 if that element of the precision matrix differs significantly with respect to the binary covariate
#' @export
#'
#'
#' @examples ls()
sig_diff_prec <- function(obj) {
  if (class(obj) != "mdine") stop("obj must be of class \"mdine\"")

  is.sig <- (obj$ci$invsigma_diff[[1]]>0) | (obj$ci$invsigma_diff[[2]]<0)

  return(1*is.sig)
}


#' Plot a network for each group
#'
#' @param obj An object of class \code{mdine}
#'
#' @return No return value
#' @export
#'
#'
#' @examples ls()
plot_networks <- function(obj) {
  if (class(obj) != "mdine") stop("obj must be of class \"mdine\"")
  return(1)
}


#' Create an \code{igraph} object from weighted adjacency matrix
#'
#' @param w.adj The weighted adjacency matrix to be converted to an \code{igraph} object
#' @param v.col The vertex colours.  Numeric vector of length equal to number of OTUs in network
#' @param e.col The edge colour.  If NULL, then colour green/red based on positive/negative weights.
#'
#' @return An object of class \code{igraph} corresponding to the provided adjacency matrix
#' @export
#'
#' @importFrom igraph categorical_pal graph.adjacency E V
#'
#' @examples  ls()
adj2ig <- function(w.adj=NULL, v.col=NULL, e.col=NULL) {
  J <- NCOL(w.adj)

  if (is.null(col)) {
    palette <- igraph::categorical_pal(J)
  } else {
    palette <- col
  }

  graph = graph.adjacency(w.adj, mode="undirected", weighted=TRUE, diag = FALSE)
  if (is.null(e.col)) {
    igraph::E(graph)$color[igraph::E(graph)$weight>0] <- "forestgreen"
    igraph::E(graph)$color[igraph::E(graph)$weight<=0] <- "orangered"
  }

  igraph::V(graph)$color <- palette

  return(graph)
}



