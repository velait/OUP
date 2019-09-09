data {
  int<lower=0> N_times;             //number of time points
  int<lower=0> N_series;            //Number of series
  real time[N_times];               //observation times
  matrix[N_series, N_times] Y;      //observations
}

transformed data {
  
  vector [N_series] mu = rep_vector(0, N_series);
  
}


parameters {
  
  vector<lower=0> [N_series] length_scale;        // mean-reversion
  vector<lower=0> [N_series] stat_var;             // volatility
  // vector [N_series] mu;                       // long-term mean
  // vector [N_series] mu_raw;                       // long-term mean
  
  real<lower=0> epsilon;
  
  // Hyperparameters
  real<lower=1> length_scale_a;
  real<lower=0> length_scale_b;
  // real<lower=0> inv_lambda_mean;
  // real<lower=0> inv_lambda_sd;
  // real<lower=0> sigma_a;
  // real<lower=0> sigma_b;
  real<lower=0> stat_var_a;
  real<lower=0> stat_var_b;
  // real mu_a;
  // real<lower=0> mu_b;
  
}

transformed parameters{
  
  vector<lower=0> [N_series] sigma = sqrt(rep_vector(2, N_series) .* stat_var ./ length_scale);           
  
  // vector [N_series] mu = 0 + mu_b*mu_raw;
  
  // vector<lower = 0> [N_series] kappa = rep_vector(0.5, N_series) .* square(sigma) .* inv_lambda;
  
  
  // real<lower=0> inv_lambda_mean = inv_lambda_a ./ inv_lambda_b;
  // real<lower=0> inv_lambda_sd = sqrt(inv_lambda_a) ./ inv_lambda_b;
  // real<lower=0> inv_lambda_a = square(inv_lambda_mean ./ inv_lambda_sd);
  // real<lower=0> inv_lambda_b = inv_lambda_mean ./ square(inv_lambda_sd);
  
}

model {
  
  
  for(s in 1:N_series) {
    
    Y[s, 1] ~ student_t(5, mu[s], sqrt(stat_var[s]));
    
    for(t in 2:N_times) {
      real dt = time[t] - time[t-1];
      real transition_mean = mu[s] - (mu[s] - Y[s, t-1])*exp(-dt/length_scale[s]);
      real transition_var = stat_var[s]*(1 - exp(-2*dt/length_scale[s]));
      
      Y[s, t] ~ student_t(5, transition_mean, sqrt(transition_var) + epsilon);
      
    }
    
  }
  
  
  
  // Priors **************
    length_scale ~ gamma(length_scale_a, length_scale_b);
  stat_var ~ gamma(stat_var_a, stat_var_b);
  // kappa ~ normal(kappa_a, kappa_b);
  // sigma  ~ gamma(sigma_a, sigma_b);
  // mu     ~ normal(mu_a, mu_b);
  // mu_raw ~ normal(0, 1);
  
  // Hyperperiors ********
    // inv_lambda_mean ~ normal(5, 1);
  // inv_lambda_sd ~ normal(5, 1);
  length_scale_a ~ gamma(2, 1);
  length_scale_b ~ gamma(2, 1);
  // sigma_a ~ normal(0, 1);
  // sigma_b ~ normal(0, 1);
  stat_var_a ~ gamma(2, 1);
  stat_var_b ~ gamma(2, 1);
  // kappa_a ~ normal(0, 1);
  // kappa_b ~ normal(0, 1);
  // mu_a ~ normal(0, 1);
  // mu_b ~ normal(0, 1);
  
  
  epsilon ~ normal(0, 1);
  
}



generated quantities{
  vector<lower=0> [N_series] lambda = 1 ./ length_scale;
  vector<lower=0> [N_series] unit_variance = stat_var .* (1 - exp(-2 ./ length_scale));
}