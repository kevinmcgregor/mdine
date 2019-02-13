
#' Get adjacency matrix for estimated networks based on credible intervals
#'
#' @param obj An object of class \code{mdine}
#' @param weighted Logical: should the adjacency matrices be weighted according to the corresponding partial correlation matrix?
#'
#' @return A list containing adjacency matrices for the two estimated precision matrices
#' @export
#'
#' @importFrom stats cov2cor
#'
#' @examples ls()
ci2adj <- function(obj, weighted=FALSE) {
  if (class(obj) != "mdine") stop("obj must be of class \"mdine\"")

  adj0 <- (obj$ci$invsigma0[[1]]>0) | (obj$ci$invsigma0[[2]]<0)
  adj1 <- (obj$ci$invsigma1[[1]]>0) | (obj$ci$invsigma1[[2]]<0)

  if (weighted) {
    adj0 <- -stats::cov2cor(obj$post_mean$invsigma0)*adj0
    adj1 <- -stats::cov2cor(obj$post_mean$invsigma1)*adj1
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


#' Plot a network (estimated from \code{mdine}) for each group
#'
#' @param obj An object of class \code{mdine}
#' @param v.col Vertex colours. If null, a colour blind-friendly pallete is used
#' @param e.col Edge colour.
#' @param lay0 \code{igraph} layout for group 0
#' @param lay1 \code{igraph} layout for group 1
#' @param lab0 Main label for group 0 network
#' @param lab1 Main label for group 1 network
#' @param scale_line_width Scaling factor for with of network edges
#' @param vertex.size Scaling factor for vertex size
#' @param vertex.labs Character vector containing vertex labels
#' @param vertex.label.cex Scaling factor for vertex labels
#'
#' @export
#'
#' @importFrom igraph plot.igraph E layout_in_circle
#' @importFrom graphics layout legend par plot.new
#'
#' @examples ls()
plot_networks <- function(obj, v.col=NULL, e.col=NULL, lay0=layout_in_circle, lay1=layout_in_circle,
                          lab0="Group 0", lab1="Group 1", scale_line_width=30, vertex.size=50,
                          vertex.labs=NULL, vertex.label.cex=NULL) {
  if (class(obj) != "mdine") stop("obj must be of class \"mdine\"")

  w.adj <- ci2adj(obj, weighted = TRUE)

  g0 <- adj2ig(w.adj$adj0, v.col, e.col)
  g1 <- adj2ig(w.adj$adj1, v.col, e.col)

  if (is.function(lay0)) {
    l0 = lay0(g0)
  } else {
    l0 = lay0
  }

  if (is.function(lay1)) {
    l1 = lay1(g1)
  } else {
    l1 = lay1
  }

  if (is.null(igraph::E(g0)$weight)) {
    lwd0 <- 0
  } else {
    lwd0 <- abs(igraph::E(g0)$weight)*scale_line_width
  }

  if (is.null(igraph::E(g1)$weight)) {
    lwd1 <- 0
  } else {
    lwd1 <- abs(igraph::E(g1)$weight)*scale_line_width
  }

  layout(matrix(c(1,2,3,3), ncol=2, byrow=TRUE), heights=c(5, 2), widths=c(6,6,1.5))
  par(mai=rep(0.2, 4))
  igraph::plot.igraph(g0, layout=l1,
       edge.width=lwd0, main=lab0,
       vertex.size=vertex.size, vertex.shape="circle",
       vertex.label=vertex.labs, vertex.label.color="black", vertex.label.cex=vertex.label.cex)
  par(mai=rep(0.2, 4))
  igraph::plot.igraph(g1, layout=l1,
       edge.width=lwd1, main=lab1,
       vertex.size=vertex.size, vertex.shape="circle",
       vertex.label=vertex.labs, vertex.label.color="black", vertex.label.cex=vertex.label.cex)
  plot.new()
  legend("center", legend=c("Positive assoc.", "Negative assoc."),
         lty=c(1,1), col=c("forestgreen","orangered"), ncol=2, cex=1.3, lwd = 6, box.col="white")

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

  pal <- igraph::categorical_pal(J)

  graph = graph.adjacency(w.adj, mode="undirected", weighted=TRUE, diag = FALSE)
  if (is.null(e.col)) {
    igraph::E(graph)$color[igraph::E(graph)$weight>0] <- "forestgreen"
    igraph::E(graph)$color[igraph::E(graph)$weight<=0] <- "orangered"
  }

  igraph::V(graph)$color <- pal

  return(graph)
}



