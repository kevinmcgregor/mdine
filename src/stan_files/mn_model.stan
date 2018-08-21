// Stan multinomial normal model with 2 covariance matrices and trying Cholesky decomp on covar matrices.  Penalty on diagonal of prec matrices.  Also using MLE for lambda as mean of exponential prior.
// Includes Frobenius norm on lower triangular part of prec matrices

functions {
  //Frobenius norm of off diagonal lower-triangular elements of matrix
  real frobenius_lower(matrix y) {
    int mat_dim;
    real f;
    f = 0;
    mat_dim = rows(y);
    for (i in 1:mat_dim) {
      for (j in 1:(i-1)) {
        f = f + fabs(y[i,j]);
      }
    }
    return f;
  }
  real nat_con(matrix y) {
    return log(mean(exp(eigenvalues_sym(fabs(y)))));
  }
}
data {
  int<lower=0> n;
  int<lower=0> k;
  int<lower=0> p;
  int<lower=0> counts[n,k];
  real<lower=0> lam_mle;
  vector[n] status;
  vector[n] offset;
  matrix[n,p] covars;
}
transformed data {
  real<lower=0> n_real;
  matrix[n,p] Q_ast;
  matrix[p,p] R_ast;
  matrix[p,p] R_ast_inverse;
  real<lower=0> inv_lam_mle;
  n_real = n;
  Q_ast = qr_Q(covars)[, 1:p] * sqrt(n_real - 1);
  R_ast = qr_R(covars)[1:p, ] / sqrt(n_real - 1);
  R_ast_inverse = inverse(R_ast);
  inv_lam_mle = 1/lam_mle;
}
parameters {
  real<lower=0> lambda;
  matrix[p,k-1] theta;
  vector[k-1] lin_pred_rand[n];
  cholesky_factor_cov[k-1] L0;
  cholesky_factor_cov[k-1] L1;
}
transformed parameters {
  //cov_matrix[k-1] sigma0;
  //cov_matrix[k-1] sigma1;
  //cov_matrix[k-1] invsigma0;
  //cov_matrix[k-1] invsigma1;
  cov_matrix[k-1] invsigma0;
  cov_matrix[k-1] invsigma1;
  simplex[k] probs[n];
  real<lower=0> invlambda;
  real<lower=0> lambda_half;
  vector[k-1] lin_pred[n];
  //cholesky_factor_cov[k-1] L0_inv;
  //cholesky_factor_cov[k-1] L1_inv;
  for (i in 1:n) {
    lin_pred[i] = (covars[i]*theta+offset[i])';
    probs[i] = softmax(append_row(lin_pred_rand[i], 0));
  }

  // sigma0 = inverse(invsigma0);
  // sigma1 = inverse(invsigma1);
  // L0 = cholesky_decompose(sigma0);
  // L1 = cholesky_decompose(sigma1);

  //L0_inv = inverse(L0);
  //L1_inv = inverse(L1);
  //invsigma0 = inverse(multiply_lower_tri_self_transpose(L0));
  //invsigma1 = inverse(multiply_lower_tri_self_transpose(L1));
  invsigma0 = multiply_lower_tri_self_transpose(L0);
  invsigma1 = multiply_lower_tri_self_transpose(L1);

  invlambda = 1/lambda;
  lambda_half = lambda/2;
}
model {
  //Likelihood
  for (i in 1:n) {
    counts[i] ~ multinomial(probs[i]);
  }

  //Putting normal dist on linear predictor
  for (i in 1:n) {
    if (status[i]==0) {
      //lin_pred_rand[i] ~ multi_normal_cholesky(lin_pred[i], L0);
      lin_pred_rand[i] ~ multi_normal_prec(lin_pred[i], invsigma0);
    } else {
      //lin_pred_rand[i] ~ multi_normal_cholesky(lin_pred[i], L1);
      lin_pred_rand[i] ~ multi_normal_prec(lin_pred[i], invsigma1);
    }
  }

  //Laplace priors for elements of precision matrices
  for (j1 in 1:(k-1)) {
    //for (j2 in 1:(j1-1)) {
      for (j2 in 1:j1) {
      if (j1==j2) {
        invsigma0[j1,j2] ~ exponential(lambda_half);
        invsigma1[j1,j2] ~ exponential(lambda_half);
      } else {
        invsigma0[j1,j2] ~ double_exponential(0,invlambda);
        invsigma1[j1,j2] ~ double_exponential(0,invlambda);
      }
    }
  }
  //Non-informative normal priors for fixed effect params
  for (j1 in 1:p) {
    for (j2 in 1:(k-1)) {
      theta[j1,j2] ~ normal(0,10000);
    }
  }

  //Prior for lambda
  lambda ~ exponential(inv_lam_mle);
}
generated quantities {
  matrix[p,k-1] beta;
  matrix[k-1,k-1] invsigma_diff;
  real frob;
  real natcon;
  //Transform fixed effects back to original scale
  beta = R_ast_inverse * theta;
  invsigma_diff = invsigma1 - invsigma0;
  frob = frobenius_lower(invsigma1)-frobenius_lower(invsigma0);
  natcon = nat_con(invsigma1)-nat_con(invsigma0);
}
