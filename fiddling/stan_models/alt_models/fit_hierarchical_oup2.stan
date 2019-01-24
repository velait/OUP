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
  
}

transformed data {
  vector[N] mu = rep_vector(0, N);
}

parameters {
  vector<lower=0> [S] inv_rho;
  vector<lower=0> [S] alpha_raw;
  real<lower=0> sigma;
  
  
  // hyperparameters
  real<lower=0> alpha_mean;
  real<lower=0> alpha_sd;
  // real<lower=0> inv_rho_shape;
  // real<lower=0> inv_rho_rate;
  
  real<lower=0> inv_rho_mean;
  real<lower=0> inv_rho_var;
  
}

transformed parameters{
  
  // non-center alpha
  vector<lower=0> [S] alpha = alpha_mean + alpha_sd*alpha_raw;
  
  real<lower=0> inv_rho_shape = square(inv_rho_mean)/inv_rho_var;
  real<lower=0> inv_rho_rate = inv_rho_mean/inv_rho_var;
  
  // real<lower=0> inv_rho_mean = inv_rho_shape/inv_rho_rate;
  // real<lower=0> inv_rho_var = inv_rho_shape/square(inv_rho_rate);
  
  vector<lower=0> [S] rho = 1 ./ inv_rho;

}
model {
  
  for(i in 1:S) {
    
    // real rho = pow(inv_rho[i], -1);
    
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
  inv_rho ~ gamma(inv_rho_shape, inv_rho_rate);
  // alpha ~ normal(alpha_mean, alpha_sd);
  alpha_raw ~ normal(0, 1);
  sigma ~ normal(0, 1);
  
  // hyperpriors
  
  // inv_rho_shape ~ normal(50, 10);
  inv_rho_mean ~ normal(10, 5);
  inv_rho_var ~ normal(0, 2);
  
  
  alpha_mean ~ normal(0, 10);
  alpha_sd ~ normal(0, 10);
  
}
