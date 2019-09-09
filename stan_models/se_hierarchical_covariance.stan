functions {
  
  // compute shape matrix 
  matrix cov_exp(real var_par, real lengthscale, real[] time, int T){
    matrix[T, T] covariance;
    
    
    for(j in 1:T) {
      for(i in 1:T) {
        covariance[i, j] = square(var_par)*exp(-fabs(time[i]-time[j])/square(lengthscale));
      }
    }
  
    return covariance;
  }
  
}


data {
  int<lower=1> N;           //number of observation
  int<lower=1> N_series;
  real x[N];                
  matrix[N, N_series] y;                //observations
}

transformed data {
  vector[N] mu = rep_vector(0, N);
}

parameters {
  
  real<lower=0> rho[N_series];        // lengthscale
  real<lower=0> alpha[N_series];      // standard deviation
  real<lower=0> sigma[N_series];      // measurement error std
  
  
  real<lower=0> rho_a;
  real<lower=0> rho_b;
  real<lower=0> alpha_a;
  real<lower=0> alpha_b;
  
  
}

transformed parameters {
  
  
  
}

model {
  
  
  for(i in 1:N_series) {
  
  matrix[N, N] L_K;
  // matrix[N, N] K = cov_exp_quad(x, alpha[i], rho[i]);
  matrix[N, N] K =  cov_exp(alpha[i], rho[i], x, N);
  real sq_sigma = square(sigma[i]);
  
  // diagonal elements
  for (n in 1:N)
    K[n, n] = K[n, n] + sq_sigma;
  
  L_K = cholesky_decompose(K);

  y[, i] ~ multi_normal_cholesky(mu, L_K);
    
  }

  
    
  rho ~ gamma(rho_a, rho_b);
  alpha ~ gamma(alpha_a, alpha_b);
  sigma ~ normal(0, 1);
  
  rho_a ~ gamma(2, .1);
  rho_b ~ gamma(2, .1);
  alpha_a ~ gamma(2, .1);
  alpha_b ~ gamma(2, .1);
  
}


generated quantities {
  
  real length_scale[N_series] = rho;   
  real stat_var[N_series] = square(alpha);
}
