functions {
  
  // compute shape matrix 
  matrix exp_cov_sq(real nu, real kappa, real inv_lambda, real[] time, int T){
    matrix[T, T] covar;
    
    
    for(i in 1:T) {
      for(j in 1:T) {
        covar[i, j] = square(kappa)*exp(-fabs(time[i]-time[j])/(2*square(inv_lambda)));
      }
    }
    
    covar = ((nu - 2)/nu)*covar;
    
    return covar;
  }
  
  // make vector with equal inputs
  vector vectorize(real mu, int T) {
    vector[T] vec;
    for(i in 1:T) {
      vec[i] = mu;
    }
    return(vec);
  }
  
}

data {
  int<lower=0> N_times;           //number of time points
  int<lower=0> N_series;           //number of series
  matrix[N_series, N_times] Y;           //observations
  real time[N_times];             //observation times
  
  // real kappa_prior[N_series, 2];
  
}

transformed data {
  
  real<lower=0> nu = 5;
  
}

parameters {
  
  vector[N_series] log_inv_lambda;        // mean-reversion
  vector[N_series] log_kappa;             // volatility
  vector [N_series] mu_raw;               // long-term mean
  
  // Hyperparameters
  real<lower=1> inv_lambda_a;
  real<lower=0> inv_lambda_b;
  
  // real<lower=0> inv_lambda_mu;
  // real<lower=0> inv_lambda_psi;
  
  
  real<lower=0> kappa_a;
  real<lower=0> kappa_b;
  
  // real<lower=0> kappa_mu;
  // real<lower=0> kappa_psi;
  
  real mu_a;
  real<lower=0> mu_b;
  
  real<lower=0> epsilon;
  
  
}

transformed parameters {
  
  vector [N_series] inv_lambda = exp(log_inv_lambda);
  vector [N_series] kappa = exp(log_kappa);
  
  // non-center mu and sigma
  vector [N_series] mu = mu_a + mu_b*mu_raw;
  
  // vector<lower = 0> [N] lambda = rep_vector(1, N)./inv_lambda;
  
  
  
  // real<lower=0> inv_lambda_b = inv(inv_lambda_psi);
  // real<lower=0> inv_lambda_a = inv_lambda_mu/inv_lambda_psi;
  
  // real<lower=0> kappa_a = inv(kappa_psi);
  // real<lower=0> kappa_b = kappa_mu/kappa_psi;
  
  // mode and variance of gamma prior
  // real<lower=0> inv_gamma_mode = lambda_sd/(lambda_mean + 1);
  // real<lower=0> inv_gamma_variance = (lambda_sd^2)/((lambda_mean-1)^2*(lambda_mean-2));
}

model {
  
  // Jacobian
  target += sum(log_inv_lambda); 
  target += sum(log_kappa);
  
  for(i in 1:N_series) {
    Y[i] ~ multi_student_t(nu,
                           rep_vector(mu[i], N_times),
                           exp_cov_sq(nu, kappa[i], inv_lambda[i], time, N_times) + diag_matrix(rep_vector(square(epsilon), N_times)));
    
  }
  
  // Priors
  inv_lambda ~ inv_gamma(inv_lambda_a, inv_lambda_b);
  kappa  ~ inv_gamma(kappa_a, kappa_b);
  // kappa  ~ normal(kappa_prior[, 1], kappa_prior[, 2]);
  // mu     ~ normal(mu_a, mu_b);
  mu_raw ~ normal(0, 1);
  
  // Hyperperiors
  // inv_lambda_mean ~ normal(5, 1);
  // inv_lambda_sd ~ normal(5, 1);
  // inv_lambda_a ~ normal(2, .1);
  // inv_lambda_b ~ normal(2, .1);
  inv_lambda_a ~ gamma(2, .1);
  inv_lambda_b ~ gamma(2, 4);
  // inv_lambda_mu ~ normal(5, 1);
  // inv_lambda_psi ~ normal(1, 1);
  // sigma_a ~ normal(0, 1);
  // sigma_b ~ normal(0, 1);
  // kappa_a ~ normal(2, .1);
  // kappa_b ~ normal(2, .1);
  // kappa_mu ~ normal(1, 1);
  // kappa_psi ~ normal(1, 1);
  kappa_a ~ gamma(2, 4);
  kappa_b ~ gamma(2, .1);
  // kappa_a ~ normal(0, 1);
  // kappa_b ~ normal(0, 1);
  mu_a ~ normal(0, 1);
  mu_b ~ normal(0, 1);
  
  
  epsilon ~ normal(0, 1);
  
  // nu ~ gamma(2, .1);
  
}

generated quantities{
  vector [N_series] length_scale = inv_lambda;
  vector [N_series] stat_var = kappa;
  vector<lower=0> [N_series] lambda = 1 ./ inv_lambda;
  vector<lower=0> [N_series] sigma = sqrt(rep_vector(2, N_series) .* kappa .* lambda);
  vector<lower=0> [N_series] unit_variance = kappa .* (1 - exp(-2 ./ inv_lambda));
}