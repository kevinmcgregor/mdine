
mn_model <- function(inData, iter, chains, mc.cores, lam.fixed=FALSE, ...) {
  mod <- ifelse(lam.fixed, stanmodels$mn_model_lam_fixed, stanmodels$mn_model)
  out <- rstan::sampling(mod, data = inData, iter=iter, chains=chains, cores=mc.cores, ...)
  return(out)
}
