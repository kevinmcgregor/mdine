
mn_model <- function(inData, iter, chains, mc.cores, lam.fixed=FALSE, ...) {
  if (lam.fixed) {
    mod <- stanmodels$mn_model_lam_fixed
  } else {
    mod <- stanmodels$mn_model
  }
  out <- rstan::sampling(mod, data = inData, iter=iter, chains=chains, cores=mc.cores, ...)
  return(out)
}
