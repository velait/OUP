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
  int<lower=1> N;
  real x[N];
  vector[N] y;
  real<lower=0> est_alpha;
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
  
  // rho ~ inv_gamma(2, 5);
  rho ~ normal(0, 10);
  alpha ~ normal(est_alpha, 1);
  sigma ~ normal(0, 1);
  
  y ~ multi_normal_cholesky(mu, L_K);
}

generated quantities {
  
  real oup_sigma = sqrt(2)*alpha/rho;
  
}