
mn_model <- function(inData, iter, chains, mc.cores, ...) {
  out <- rstan::sampling(stanmodels$mn_model, data = inData, iter=iter, chains=chains, cores=mc.cores, ...)
  return(out)
}
