functions {
  
}

data {
  int<lower = 1> N;
  real x[N];  
  int<lower =1> D;
  matrix [N, D] y;
}

parameters {
  
  // parameters
  vector [D] mu;
  
  corr_matrix[D] Omega_G;
  // corr_matrix[2] Omega_B;
  vector<lower=0>[D] sigma_G;
  // vector<lower=0>[2] sigma_B;
  
  cov_matrix [D] B;
  
}

transformed parameters {
 
  cov_matrix [D] G = quad_form_diag(Omega_G, sigma_G);
  // cov_matrix [2] B = quad_form_diag(Omega_B, sigma_B);
  
}


model {
  
  y[1, ] ~ multi_normal(mu, G);
  
  
  for (i in 2:N) {

  real delta_t = x[i] - x[i-1];


  vector [D] m = mu + matrix_exp(-B*delta_t) * (to_vector(y[i-1, ]) - mu);
  matrix [D, D] v = G - matrix_exp(-B*delta_t) * G * matrix_exp(-B*delta_t);

  y[i, ] ~ multi_normal(m, v);

  }
  

  mu ~ normal(0, 2);
  
  sigma_G ~ cauchy(0, 5);
  // sigma_B ~ cauchy(0, 5);
  
  Omega_G ~ lkj_corr(1);
  // Omega_B ~ lkj_corr(1);

  // G[1, 1] ~ normal(0, 1);
  // G[2, 2] ~ normal(0, 1);
  // G[1, 2] ~ normal(0, 1);
  // G[2, 1] ~ normal(0, 1);
  
  for(i in 1:D) {
    for(j in 1:D) {
      if(i == j) {
        B[i, i] ~ gamma(4, 10);
      } else {
        B[i, j] ~ normal(0, 1);
      }
    }
  }

}