
extract_post_mean <- function(stan.fit, n, k, p, lam.null) {
  post_mean <- NULL
  post_mean$beta <- matrix(rstan::get_posterior_mean(stan.fit, "beta")[,"mean-all chains"], p, k-1, byrow=TRUE)
  post_mean$invsigma0 <- matrix(rstan::get_posterior_mean(stan.fit, "invsigma0")[,"mean-all chains"], k-1, k-1, byrow=TRUE)
  post_mean$invsigma1 <- matrix(rstan::get_posterior_mean(stan.fit, "invsigma1")[,"mean-all chains"], k-1, k-1, byrow=TRUE)
  post_mean$invsigma_diff <- matrix(rstan::get_posterior_mean(stan.fit, "invsigma_diff")[,"mean-all chains"], k-1, k-1, byrow=TRUE)
  post_mean$w <- matrix(rstan::get_posterior_mean(stan.fit, "lin_pred_rand")[,"mean-all chains"], n, k-1, byrow=TRUE)
  post_mean$probs <- matrix(rstan::get_posterior_mean(stan.fit, "probs")[,"mean-all chains"], n, k-1, byrow=TRUE)
  post_mean$frob <- rstan::get_posterior_mean(stan.fit, "frob")[,"mean-all chains"]
  post_mean$natcon <- rstan::get_posterior_mean(stan.fit, "natcon")[,"mean-all chains"]
  if (lam.null) {
    post_mean$lambda <- rstan::get_posterior_mean(stan.fit, "lambda")[,"mean-all chains"]
  } else {
    get_posterior_mean$lambda <- NULL
  }

  return(post_mean)
}

extract_ci <- function(stan.fit,n,k,p,lam.null,quant) {
  col <- paste0(as.character(quant*100),"%")
  param_smry <- rstan::summary(stan.fit, pars=c("beta","lin_pred_rand","invsigma0","invsigma1","invsigma_diff",
                                        "frob","natcon"))$summary
  beta_lower <- matrix(param_smry[grep("beta",rownames(param_smry)),col[1]],p,k-1,byrow=TRUE)
  beta_upper <- matrix(param_smry[grep("beta",rownames(param_smry)),col[2]],p,k-1,byrow=TRUE)
  w_lower <- matrix(param_smry[grep("lin_pred_rand",rownames(param_smry)),col[1]],n,k-1,byrow=TRUE)
  w_upper <- matrix(param_smry[grep("lin_pred_rand",rownames(param_smry)),col[2]],n,k-1,byrow=TRUE)
  probs_lower <- matrix(param_smry[grep("probs",rownames(param_smry)),col[1]],n,k,byrow=TRUE)
  probs_upper <- matrix(param_smry[grep("probs",rownames(param_smry)),col[2]],n,k,byrow=TRUE)
  invsigma0_lower <- matrix(param_smry[grep("invsigma0",rownames(param_smry)),col[1]],k-1,k-1,byrow=TRUE)
  invsigma0_upper <- matrix(param_smry[grep("invsigma0",rownames(param_smry)),col[2]],k-1,k-1,byrow=TRUE)
  invsigma1_lower <- matrix(param_smry[grep("invsigma1",rownames(param_smry)),col[1]],k-1,k-1,byrow=TRUE)
  invsigma1_upper <- matrix(param_smry[grep("invsigma1",rownames(param_smry)),col[2]],k-1,k-1,byrow=TRUE)
  invsigma_diff_lower <- matrix(param_smry[grep("invsigma_diff",rownames(param_smry)),col[1]],k-1,k-1,byrow=TRUE)
  invsigma_diff_upper <- matrix(param_smry[grep("invsigma_diff",rownames(param_smry)),col[2]],k-1,k-1,byrow=TRUE)
  frob_lower <- param_smry[grep("frob",rownames(param_smry)),col[1]]
  frob_upper <- param_smry[grep("frob",rownames(param_smry)),col[2]]
  natcon_lower <- param_smry[grep("natcon",rownames(param_smry)),col[1]]
  natcon_upper <- param_smry[grep("natcon",rownames(param_smry)),col[2]]
  if (lam.null) {
    lambda_lower <- rstan::summary(stan.fit, pars=c("lambda"))$summary[,col[1]]
    lambda_upper <- rstan::summary(stan.fit, pars=c("lambda"))$summary[,col[2]]
  } else {
    lambda_lower <- NULL
    lambda_upper <- NULL
  }

  return(list(
    beta = list(lower=beta_lower, upper=beta_upper),
    invsigma0 = list(lower=invsigma0_lower, upper=invsigma0_upper),
    invsigma1 = list(lower=invsigma1_lower, upper=invsigma1_upper),
    invsigma_diff = list(lower=invsigma_diff_lower, upper=invsigma_diff_upper),
    lambda = c(lambda_lower, lambda_upper),
    frob = c(frob_lower, frob_upper),
    natcon = c(natcon_lower, natcon_upper) ))
}

