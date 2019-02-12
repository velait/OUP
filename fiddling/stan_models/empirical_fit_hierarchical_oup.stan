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
  
  real<lower=0> est_alpha_mean;
  real<lower=0> est_alpha_sd;
  
}
transformed data {
  vector[N] mu = rep_vector(0, N);
}
parameters {
  vector<lower=0> [S] rho;
  vector<lower=0> [S] alpha;
  real<lower=0> sigma;
  
  
  // hyperparameters
  real<lower=0> rho_shape;
  real<lower=0> rho_rate;
  // real<lower=0> alpha_sd;
  
  
}

transformed parameters{
  
  // vector<lower=0> [S] alpha = est_alpha + sqrt(est_alpha_var)*alpha_raw;
  
}
model {
  
  for(i in 1:S) {
    matrix[N, N] L_K;
    matrix[N, N] K = cov_exp_abs(x, alpha[i], rho[i], N);
    real sq_sigma = square(sigma);
    
    // diagonal elements
    for (n in 1:N) {
      K[n, n] = K[n, n] + sq_sigma;
    }
    
    
    L_K = cholesky_decompose(K);
    
    y[i] ~ multi_normal_cholesky(mu, L_K);
  }
  
  
  // priors
  rho ~ inv_gamma(rho_shape, rho_rate);
  alpha ~ normal(est_alpha_mean, est_alpha_sd);
  // alpha_raw ~ normal(0, 1);
  sigma ~ normal(0, 1);
  
  // hyperpriors
  rho_shape ~ normal(10, 5);
  rho_rate ~ normal(50, 5);
  // alpha_mean ~ normal(est_alpha, est_alpha_var);
  // alpha_sd ~ normal(0, 1);
  // 
}
