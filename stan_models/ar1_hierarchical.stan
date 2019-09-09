data {
  int<lower=0> N_times;             //number of time points
  int<lower=0> N_series;            //Number of series
  real time[N_times];               //observation times
  matrix[N_series, N_times] Y;      //observations
}

parameters {
  
  // vector<lower=0, upper = 1> [N_series] phi;        // mean-reversion
  vector<lower=0> [N_series] sigma;      // volatility
  // vector [N_series] c;                   // long-term mean
  
  vector<lower=0, upper = 1> [N_series] lambda;
  vector [N_series] mu_raw;
  // vector [N_series] sigma_raw;
  
  // Hyperparameters
  real<lower=1> lambda_a;
  real<lower=0> lambda_b;
  real<lower=0> sigma_a;
  real<lower=0> sigma_b;
  real mu_a;
  real<lower=0> mu_b;
  
}

transformed parameters{

  vector [N_series] mu = mu_a + mu_b*mu_raw;
  // vector [N_series] sigma = sigma_a + sigma_b*sigma_raw;
  
  vector [N_series] stat_mu = (lambda .* mu) ./ (rep_vector(1.0, N_series) + lambda);
  vector<lower=0> [N_series] stat_variance = square(sigma) ./ (1 - square(lambda));
}

model {
  
  
  for(s in 1:N_series) {
    
      Y[s, 1] ~ normal(stat_mu[s], sqrt(stat_variance[s]));
  
  for(t in 2:N_times) {
    
    Y[s, t] ~ normal(lambda[s]*mu[s] - lambda[s]*Y[s, t-1], square(sigma[s]));
    
  }
    
    
  }

  
  // Priors
  lambda ~ gamma(lambda_a, lambda_b);
  sigma  ~ gamma(sigma_a, sigma_b);
  // mu     ~ normal(mu_a, mu_b);
  mu_raw ~ normal(0, 1);
  // sigma_raw ~ normal(0, 1);
  
  // Hyperpriors
  lambda_a ~ normal(5, 2);
  lambda_b ~ normal(5, 2);
  sigma_a ~ normal(5, 2);
  sigma_b ~ normal(5, 2);
  // sigma_a ~ normal(0, 1);
  // sigma_b ~ normal(0, 1);
  mu_a ~ normal(0, 1);
  mu_b ~ normal(0, 1);
  
}

generated quantities {
  
  vector [N_series] phi = -lambda;
  vector [N_series] c = lambda .* mu;
}



