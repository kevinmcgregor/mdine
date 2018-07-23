#' Microbiome differential network estimation
#'
#' @export
#'
#' @param ... Other arguments passed to `rstan::sampling`
#' @param Y The (unrarefied) taxa count matrix with rows as samples and columns as taxa.  The last column is
#' the reference category, and is not included in the estimated network.
#' @param X The model matrix (including an intercept column)
#' @param Z Binary variable over which the network is assumed to vary.  If NULL then a single network
#' is assumed.
#' @param lambda Network penalization parameter
#' @param offset A vector containing an offset term for each subject
#' @param mc.cores The number of cores to run MCMC chains in parallel
#' @param iter The number of MCMC iterations.  By default the first half of the iterations will be used as warmup.
#' @param chains The number of MCMC chains.
#' @param quant Lower and upper quantiles of the posterior distribution to create credible intervals.
#' @param method Whether to solve using MCMC sampling, or through maximum a posteriori estimation (MAP)
#'
#' @return A list containing posterior means for the model parameters, credible intervals,
#' and the stanfit object.
#'
#' @examples ls()
mdine <- function(Y, X, Z=NULL, lambda=NULL, offset=NULL, mc.cores=1,iter=1000,
                  chains=4, quant=c(0.025, 0.975), method=c("MCMC", "MAP"), ...) {

  if (!is.matrix(Y) | any(floor(Y)!=Y) | any(Y<0)) stop("Y must be a numeric matrix containing positive counts")
  if (!is.matrix(X) | !is.numeric(X)) stop("X must be a numeric matrix")
  if (!is.vector(Z) | !is.numeric(Z) | any(!Z%in%c(0,1))) stop("Z must be a binary vector")
  if (NROW(Y)!=NROW(X) | length(Z)!=NROW(Y)) stop("Dimension mismatch in one or more of (Y,X,Z)")
  if (!is.null(lambda)) {if (lambda<=0 | !is.numeric(lambda)) stop("lambda must be positive numeric")}

  n <- NROW(Y)
  k <- NCOL(Y)
  p <- NCOL(X)

  method <- match.arg(method)

  if (method=="MCMC") {
    lam_mle <- get_lam_mle(Y, X, Z)
    if (is.null(offset)) offset <- rep(0, n)
    s.data <- list(counts=Y, covars=X, status=Z, n=n, k=k, p=p, lam_mle=lam_mle,
                   offset=offset)
    fit <- mn_model(s.data, iter, chains, mc.cores, ...)
    post_mean <- extract_post_mean(fit, is.null(lambda))
    ci <- extract_ci(fit, is.null(lambda), quant=quant)
  } else if (method=="MAP") {
    #TODO
  }

  if (!is.null(lambda)) post_mean$lambda <- lambda

  return(list(stan.fit=fit, post_mean=post_mean, ci=ci))
}

