functions {
  
  matrix gen_oup(real[] x, vector sigma, vector rho){
    
    int T = size(x);
    
    matrix[T, T] covariance;
    
    for(i in 1:T) {
      for(j in 1:T) {
        real diff_sq = square(fabs(x[i]-x[j]));
        real mean_rho = (rho[i] + rho[j])/2;
        real sqrt_Q = diff_sq/sqrt(mean_rho);
        
        covariance[i, j] = sigma[i]*sigma[j]*pow(rho[i]*rho[j], .25)*exp(-sqrt_Q)/sqrt(mean_rho);
        
        
        
        
        
        // real diff_sq = square(fabs(x[i]-x[j]));
        // real rho_sq_plus = square(rho[i]) + square(rho[j]);
        // 
        // covariance[i, j] = sigma[i]*sigma[j];
        // 
        // covariance[i, j] *= sqrt((2*rho[i]*rho[j])/(rho_sq_plus));
        // 
        // covariance[i, j] *= exp(-diff_sq/rho_sq_plus);
        
        
      }
    }
    
    return covariance;
  }
  
}

data {
  int<lower=1> N;
  real x[N];
  vector[N] y; // Observations
  
  // int<lower=1> N_pred; // prediction grid resolution
  // real x_pred[N_pred]; // Prediction time points
  
  // Parameter process hyperparameters -- note, that these are given, not learned
  real<lower=0> alpha_rho;
  real<lower=0> beta_rho;
  real<lower=0> mu_rho;
  
  real<lower=0> alpha_sigma;
  real<lower=0> beta_sigma;
  real<lower=0> mu_sigma;
  
  
}

transformed data {
  
  # Process means
  vector[N] zeros;
  vector[N] log_rho_zeros;
  vector[N] log_sigma_zeros;
  
  log_rho_zeros = rep_vector(mu_rho, N);
  log_sigma_zeros = rep_vector(mu_sigma, N);
  zeros = rep_vector(0, N);
  
}

parameters {
  // Processes latent vectors (whitening tms.)
  vector[N] eta_f;
  vector[N] eta_log_rho;
  vector[N] eta_log_sigma;
  
  real<lower=0> epsilon;
}

transformed parameters {
  
  vector[N] f;
  vector[N] rho;
  vector[N] log_rho;
  vector[N] sigma;
  vector[N] log_sigma;
  
  {
    
    // rho
    matrix[N, N] L_log_rho;
    matrix[N, N] K_log_rho;
    
    // sigma
    matrix[N, N] L_log_sigma;
    matrix[N, N] K_log_sigma;
    
    // f
    matrix[N, N] L;
    matrix[N, N] K;
    
    
    
    K_log_rho = cov_exp_quad(x, alpha_rho, beta_rho);
    for (n in 1:N)
      K_log_rho[n, n] = K_log_rho[n, n] + 1e-12;
    L_log_rho = cholesky_decompose(K_log_rho);
    log_rho = L_log_rho * eta_log_rho;
    rho = exp(log_rho);
    
    
    K_log_sigma = cov_exp_quad(x, alpha_sigma, beta_sigma);
    for (n in 1:N)
      K_log_sigma[n, n] = K_log_sigma[n, n] + 1e-12;
    L_log_sigma = cholesky_decompose(K_log_sigma);
    log_sigma = L_log_sigma * eta_log_sigma;
    sigma = exp(log_sigma);
    
    
    
    K = gen_oup(x, sigma, rho);
    for (n in 1:N)
      K[n, n] = K[n, n] + 1e-12;
    L = cholesky_decompose(K);
    f = L * eta_f;
  }
  
}

model {
  
  
  eta_log_rho ~ normal(0, 2);
  eta_log_sigma ~ normal(0, 2);
  eta_f ~ normal(0, 2);
  epsilon ~ normal(0, 1);
  
  y ~ normal(f, epsilon);
  
}

