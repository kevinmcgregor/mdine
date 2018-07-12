#' Run multinomial normal model through STAN
#'
#'
#' @description A function to run the multinomial normal model through STAN.  This is not meant to be called
#' directly (use `mdine` instead).
#' @param inData List containing all necessary elements for the input to the model
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#'
#' @return An object of class `stanfit` containing MCMC samples for the multinomial normal model
#'
mn_model <- function(inData, iter, chains, mc.cores, ...) {
  out <- rstan::sampling(stanmodels$mn_model, data = inData, iter=iter, chains=chains, cores=mc.cores, ...)
  return(out)
}
