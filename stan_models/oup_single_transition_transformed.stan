data {
  int<lower=0> T;             //number of time points
  vector[T] time;               //observation times
  vector[T] Y;                //observations
}

transformed data {
  
  real max_dt = time[T] - time[1];
  real mean_Y = mean(Y);
  real sd_Y = sd(Y);
  
  
  // Scale observations and times
  vector[T] sc_time = time ./ max_dt;
  vector[T] sc_T = (Y - mean_Y) ./ sd_Y;
  
}

parameters {
  
  real<lower=0> sc_inv_lambda;        // mean-reversion
  real<lower=0> sc_kappa;           // volatility
  real sc_mu;                       // long-term mean
  
  // real epsilon;
  
}

transformed parameters{
  
  real<lower = 0> sc_sigma = sqrt(2*sc_kappa/sc_inv_lambda);
  
}

model {
  
  Y[1] ~ normal(sc_mu, sc_kappa);
  
  for(i in 2:T) {
    
    real transition_mean = sc_mu - (sc_mu - Y[i-1])*exp(-pow(sc_inv_lambda, -1));
    real transition_var = sc_kappa*(1 - exp(-2*pow(sc_inv_lambda, -1)));
    
    Y[i] ~ normal(transition_mean, sqrt(transition_var));
    
  }
  
  // Priors
  sc_inv_lambda ~ normal(0, 1);
  sc_kappa ~ normal(0, 1);
  
  
  
  sc_inv_lambda ~ gamma(2, .5*max_dt);
  sc_kappa  ~ gamma(2, 2*sd_Y);
  sc_mu     ~ normal(0, 1);
  
}



generated quantities{
  real sc_lambda = pow(sc_inv_lambda, -1);
  real sc_unit_variance = sc_kappa*(1 - exp(-2*sc_inv_lambda));
  
  
  // Unscale the parameters
  real lambda = sc_lambda*max_dt;
  real mu = sc_mu;
  real kappa = sc_kappa*sd_Y;
  real sigma = sqrt(kappa*lambda);
  
}