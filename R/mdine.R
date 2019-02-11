#' Microbiome differential network estimation
#'
#' @export
#' @aliases mdine
#' @param ... Other arguments passed to `rstan::sampling`
#' @param Y The (unrarefied) taxa count matrix with rows as samples and columns as taxa.  The last column is
#' the reference category, and is not included in the estimated network.
#' @param X The model matrix (including an intercept column)
#' @param Z Binary variable over which the network is assumed to vary.
#' @param lambda Network penalization parameter.  If NULL, then lambda is estimated
#' @param offset A vector containing an offset term for each subject
#' @param mc.cores The number of cores to run MCMC chains in parallel
#' @param iter The number of MCMC iterations.  By default the first half of the iterations will be used as warmup.
#' @param chains The number of MCMC chains.
#' @param quant Lower and upper quantiles of the posterior distribution to create credible intervals.
#'
#' @return A an object of class "mdine" containing posterior means for the model parameters, credible intervals,
#' and the stanfit object.
#'
#' @examples ls()
mdine <- function(Y, X, Z, lambda=NULL, offset=NULL, mc.cores=1,iter=1000,
                  chains=4, quant=c(0.025, 0.975), ...) {

  if (!is.matrix(Y) | any(floor(Y)!=Y) | any(Y<0)) stop("Y must be a numeric matrix containing positive counts")
  if (!is.matrix(X) | !is.numeric(X)) stop("X must be a numeric matrix")
  if (!is.vector(Z) | !is.numeric(Z) | any(!Z%in%c(0,1))) stop("Z must be a binary vector")
  if (NROW(Y)!=NROW(X) | length(Z)!=NROW(Y)) stop("Dimension mismatch in one or more of (Y,X,Z)")
  if (!is.null(lambda)) {if (lambda<=0 | !is.numeric(lambda)) stop("lambda must be positive numeric")}
  if (!is.numeric(quant) | length(quant)!=2) stop("quant must be a vector of length 2")

  n <- NROW(Y)
  k <- NCOL(Y)
  p <- NCOL(X)

  if (is.null(offset)) offset <- rep(0, n)

  if (is.null(lambda)) {
    lam_mle <- get_lam_mle(Y, X, Z)
    s.data <- list(counts=Y, covars=X, status=Z, n=n, k=k, p=p, lam_mle=lam_mle,
                    offset=offset)
  } else {
    s.data <- list(counts=Y, covars=X, status=Z, n=n, k=k, p=p, lambda=lambda, offset=offset)
  }

  fit <- mn_model(s.data, iter, chains, mc.cores, !is.null(lambda), ...)

  post_mean <- extract_post_mean(fit, is.null(lambda))
  ci <- extract_ci(fit, is.null(lambda), quant=quant)

  ret <- list(stan.fit=fit, post_mean=post_mean, ci=ci)

  if (is.null(lambda)) {
    ret$lam_mle <- lam_mle
  } else {
    ret$lambda.fixed <- lambda
  }

  class(ret) <- "mdine"
  return(ret)
}

