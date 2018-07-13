
extract_post_mean = function(stan.fit, n, k, p, lam.null) {
  post_mean <- NULL
  post_mean$beta <- matrix(get_posterior_mean(stan.fit, "beta")[,"mean-all chains"], p, k-1, byrow=TRUE)
  post_mean$invsigma0 <- matrix(get_posterior_mean(stan.fit, "invsigma0")[,"mean-all chains"], k-1, k-1, byrow=TRUE)
  post_mean$invsigma1 <- matrix(get_posterior_mean(stan.fit, "invsigma1")[,"mean-all chains"], k-1, k-1, byrow=TRUE)
  post_mean$invsigma_diff <- matrix(get_posterior_mean(stan.fit, "invsigma_diff")[,"mean-all chains"], k-1, k-1, byrow=TRUE)
  post_mean$w <- matrix(get_posterior_mean(stan.fit, "lin_pred_rand")[,"mean-all chains"], n, k-1, byrow=TRUE)
  post_mean$probs <- matrix(get_posterior_mean(stan.fit, "probs")[,"mean-all chains"], n, k-1, byrow=TRUE)
  post_mean$frob <- get_posterior_mean(stan.fit, "frob")[,"mean-all chains"]
  post_mean$natcon <- get_posterior_mean(stan.fit, "natcon")[,"mean-all chains"]
  if (lam.null) {
    post_mean$lambda <- get_posterior_mean(stan.fit, "lambda")[,"mean-all chains"]
  } else {
    get_posterior_mean$lambda <- NULL
  }

  return(post_mean)
}


