
#Extract posterior means for all parameters... returns list
extract_post_mean <- function(stan.fit, lam.null) {
  pars <- c("beta", "invsigma0", "invsigma1", "invsigma_diff", "lin_pred_rand", "probs")

  post_mean <- vector("list", length(pars))
  names(post_mean) <- pars

  post_mean$lambda <- NULL
  if (lam.null) {
    pars <- c(pars, "lambda")
  }

  for (par in pars) {
    post_mean[[par]] <- get_post_mean_matrix(stan.fit, par)
  }

  return(post_mean)
}

#Extract credible intervals for all parameters... returns list
extract_ci <- function(stan.fit,lam.null,quant) {
  pars <- c("beta", "invsigma0", "invsigma1", "invsigma_diff", "lin_pred_rand", "probs")

  ci <- vector("list", length(pars))
  names(ci) <- pars

  ci$lambda <- NULL
  if (lam.null) {
    pars <- c(pars, "lambda")
  }

  for (par in pars) {
    ci[[par]] <- get_ci_stanfit(stan.fit, par, quant)
  }

  return(ci)
}

# Get the posterior mean of parameter and put it in matrix with appropriate dimension
#' @importFrom rstan get_posterior_mean
get_post_mean_matrix <- function(stan.fit, par) {
  dim <- stan.fit@par_dims[[par]]
  if (is.na(dim[1])) {
    param_mean <- rstan::get_posterior_mean(stan.fit, par)[,"mean-all chains"]
  } else {
    param_mean <- matrix(rstan::get_posterior_mean(stan.fit, par)[,"mean-all chains"],
                       dim[1], dim[2], byrow=TRUE)
  }
  return(param_mean)
}


# Get the credible interval of a single parameter from a stanfit and put into matrix (if appropriate)
#' @importFrom rstan summary
get_ci_stanfit <- function(stan.fit, par, quant) {
  col <- 4:(3+length(quant))
  param_smry <- rstan::summary(stan.fit, pars=par, probs=quant)$summary[,col]
  dim <- stan.fit@par_dims[[par]]

  if (!is.na(dim[1])) {
    ret <- vector("list", length(quant))
    names(ret) <- paste0(100*quant, "%")
    for (q in 1:length(quant)) {
      ret[[q]] <- matrix(param_smry[,q], dim[1], dim[2])
    }
  } else {
    ret <- param_smry
  }

  return(ret)
}

