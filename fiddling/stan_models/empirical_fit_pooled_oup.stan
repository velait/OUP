// EMPIRICAL POOLED FITTER

functions {
  matrix cov_exp_abs(real[] time, real alpha, real rho, int N){
    matrix[N, N] covariance;
    
    
    for(i in 1:N) {
      for(j in 1:N) {
        covariance[i, j] = square(alpha)*exp(-fabs(time[i]-time[j])/square(rho));
      }
    }
    
    return covariance;
  }
}

data {
  int<lower=1> N; // number of points
  int<lower=1> S; // number of series
  real x[N];      // time
  matrix[S, N] y;    // observations
  
  real<lower=0> est_alpha_mean; // Marginal variance of data --> empirical Bayes
  real<lower=0> est_alpha_sd;
  
}
transformed data {
  vector[N] mu = rep_vector(0, N);
}
parameters {
  real<lower=0> rho;
  real<lower=0> alpha;
  real<lower=0> sigma;
}
model {
  matrix[N, N] L_K;
  matrix[N, N] K = cov_exp_abs(x, alpha, rho, N);
  real sq_sigma = square(sigma);
  
  // diagonal elements
  for (n in 1:N) {
    K[n, n] = K[n, n] + sq_sigma;
  }
  
  
  L_K = cholesky_decompose(K);
  
  rho ~ inv_gamma(4, 10);
  alpha ~ normal(est_alpha_mean, est_alpha_sd);
  sigma ~ normal(0, 1);
  
  for(i in 1:S) {
    y[i] ~ multi_normal_cholesky(mu, L_K);
  }
  
}


