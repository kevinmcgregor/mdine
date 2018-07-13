#' @import nnet
#' @import MASS
get_lam_mle <- function(counts, covar, status) {
  ref <- NCOL(counts)
  n.spec <- ref-1

  r.counts <- cbind(counts[,-ref], counts[,ref])
  fit <- nnet::multinom(r.counts~covar[,-1], trace=FALSE)
  resid <- log((counts[,-ref]+1)/(counts[,ref]+1))-tcrossprod(covar,coef(summary(fit)))

  #Covariance of residuals
  cov0 <- cov(resid[status==0,])
  cov1 <- cov(resid[status==1,])

  prec0 <- MASS::ginv(cov(resid[status==0,]))
  prec1 <- MASS::ginv(cov(resid[status==1,]))

  lt0 <- sum(abs(prec0[lower.tri(prec0)]))
  lt1 <- sum(abs(prec1[lower.tri(prec1)]))
  lam_mle <- n.spec*(2*n.spec-1)/(2*(lt0+lt1)+sum(diag(prec0))+sum(diag(prec1)))

  return(lam_mle)
}
