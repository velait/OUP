functions {
  
}

data {
  int<lower = 1> N;
  real x[N];  
  matrix [N, 2] y;
}

parameters {
  
  // parameters
  vector [2] mu1;
  vector [2] mu2;
  
  corr_matrix[2] Omega_G1;
  corr_matrix[2] Omega_G2;
  // corr_matrix[2] Omega_B;
  vector<lower=0>[2] sigma_G1;
  vector<lower=0>[2] sigma_G2;
  // vector<lower=0>[2] sigma_B;
  
  cov_matrix [2] B1;
  cov_matrix [2] B2;
  
    // latent process; mixture weight
  vector<lower = 0, upper = 1> [N] mix_weight;
  
}

transformed parameters {
  
  cov_matrix [2] G1 = quad_form_diag(Omega_G1, sigma_G1);
  cov_matrix [2] G2 = quad_form_diag(Omega_G2, sigma_G2);
  // cov_matrix [2] B = quad_form_diag(Omega_B, sigma_B);
  
}


model {
  
  target += log_mix(mix_weight[1],
              multi_normal_lpdf(y[1, ] | mu1, G1),
              multi_normal_lpdf(y[1, ] | mu2, G2));
  
  
  for (i in 2:N) {
    
    real delta_t = x[i] - x[i-1];
    
    
    vector [2] m1 = mu1 + matrix_exp(-B1*delta_t) * (to_vector(y[i-1, ]) - mu1);
    matrix [2, 2] v1 = G1 - matrix_exp(-B1*delta_t) * G1 * matrix_exp(-B1*delta_t);
    vector [2] m2 = mu2 + matrix_exp(-B2*delta_t) * (to_vector(y[i-1, ]) - mu2);
    matrix [2, 2] v2 = G2 - matrix_exp(-B2*delta_t) * G2 * matrix_exp(-B2*delta_t);
    
    target += log_mix(mix_weight[i],
              multi_normal_lpdf(y[i, ] | m1, v1),
              multi_normal_lpdf(y[i, ] | m2, v2));
              
      // Random walk prior for the latent process
      mix_weight[i] ~ normal(mix_weight[i-1], 0.25*delta_t);
      target += log_mix(mix_weight[i-1], normal_lpdf(mix_weight[i] | 1 , .25), normal_lpdf(mix_weight[i] | 0 , .25));

    
  }
  
  
  mu1 ~ normal(0, 2);
  mu2 ~ normal(0, 2);
  
  sigma_G1 ~ cauchy(0, 5);
  sigma_G2 ~ cauchy(0, 5);
  // sigma_B ~ cauchy(0, 5);
  
  Omega_G1 ~ lkj_corr(1);
  Omega_G2 ~ lkj_corr(1);
  // Omega_B ~ lkj_corr(1);
  
  // G[1, 1] ~ normal(0, 1);
  // G[2, 2] ~ normal(0, 1);
  // G[1, 2] ~ normal(0, 1);
  // G[2, 1] ~ normal(0, 1);
  // 
  B1[1, 1] ~ gamma(4, 10);
  B1[2, 2] ~ gamma(4, 10);
  B1[1, 2] ~ normal(0, 1);
  B1[2, 1] ~ normal(0, 1);
  B2[1, 1] ~ gamma(4, 10);
  B2[2, 2] ~ gamma(4, 10);
  B2[1, 2] ~ normal(0, 1);
  B2[2, 1] ~ normal(0, 1);
  
  
  
}