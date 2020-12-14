functions {
  
  matrix gen_cov_exp_quad(real[] x, vector sigma, vector rho){
    
    int T = size(x);
    
    matrix[T, T] covariance;
    
    for(i in 1:T) {
      for(j in 1:T) {
        real diff_sq = square(fabs(x[i]-x[j]));
        real rho_sq_plus = square(rho[i]) + square(rho[j]);
        
        covariance[i, j] = sigma[i]*sigma[j];
        
        covariance[i, j] *= sqrt((2*rho[i]*rho[j])/(rho_sq_plus));
        
        covariance[i, j] *= exp(-diff_sq/rho_sq_plus);
        
      }
    }
    
    return covariance;
  }
  
}

data {
  int<lower=1> N;
  real x[N];
  
  // real<lower =0> alpha_rho;
  // real<lower =0> beta_rho;
  // real<lower =0> mu_rho;
  // 
  // real<lower =0> alpha_sigma;
  // real<lower =0> beta_sigma;
  // real<lower =0> mu_sigma;
  
  vector[N] sigma;
  vector[N] rho;
  real<lower=0> epsilon;
}
transformed data {
  vector[N] zeros;
  // vector[N] log_rho_zeros;
  // vector[N] log_sigma_zeros;
  
  // log_rho_zeros = rep_vector(mu_rho, N);
  // log_sigma_zeros = rep_vector(mu_sigma, N);
  zeros = rep_vector(0, N);
  
  
  
}
model {}
generated quantities {
  
  // Generate latent rho process
  // vector[N] log_rho;
  // vector[N] rho;
  
  // vector[N] log_sigma;
  // vector[N] sigma;
  
  vector[N] y;
  vector[N] f;
  
  {
    
    // Generate lengthscale process
    // matrix[N, N] log_rho_cov;
    // matrix[N, N] L_log_rho_cov;
    // 
    // matrix[N, N] log_sigma_cov;
    // matrix[N, N] L_log_sigma_cov;
    
    matrix[N, N] cov;
    matrix[N, N] L_cov;
    
    
    // log_rho_cov = cov_exp_quad(x, alpha_rho, beta_rho);
    // 
    // for (n in 1:N)
    //   log_rho_cov[n, n] = log_rho_cov[n, n] + 1e-12;
    // 
    // L_log_rho_cov = cholesky_decompose(log_rho_cov);
    // log_rho = multi_normal_cholesky_rng(log_rho_zeros, L_log_rho_cov);
    // rho = exp(log_rho);
    // 
    
    
    // Generate variance process
    
    
    // log_sigma_cov = cov_exp_quad(x, alpha_sigma, beta_sigma);
    // 
    // for (n in 1:N)
    //   log_sigma_cov[n, n] = log_sigma_cov[n, n] + 1e-12;
    // 
    // L_log_sigma_cov = cholesky_decompose(log_sigma_cov);
    // log_sigma = multi_normal_cholesky_rng(log_sigma_zeros, L_log_sigma_cov);
    // sigma = exp(log_sigma);
    
    
    
    
    // Generate process itself
    
    
    
    
    cov = gen_cov_exp_quad(x, sigma, rho);
    
    for (n in 1:N)
      cov[n, n] = cov[n, n] + 1e-12;
    
    L_cov = cholesky_decompose(cov);
    f = multi_normal_cholesky_rng(zeros, L_cov);
  }
  for (n in 1:N)
    y[n] = normal_rng(f[n], epsilon);
  
  
}