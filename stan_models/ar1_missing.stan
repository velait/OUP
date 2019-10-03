data {
  int<lower=0> N_obs;
  int<lower=0> N_mis;
  
  int ii_obs[N_obs];
  int ii_mis[N_mis];
  real Y_obs[N_obs];
}

transformed data {
  
  // real epsilon = 0.1; // error standard deviation
  int N = N_obs + N_mis;
  
}

parameters {
  
  real mu;
  real<lower=0> sigma;           // volatility
  real<lower=-1, upper=1> lambda;        // mean-reversion
  real<lower = 0> epsilon;
  
  
  // real Y_mis[N_mis];
  
  vector[N_obs + N_mis] latent_Y;
}

transformed parameters{
  real stat_mu = (lambda*mu)/(1 + lambda);
  real<lower=0> stat_var = (sigma^2)/(1 - lambda^2);
  
  
  // Y[ii_obs] = Y_obs;
  // Y[ii_mis] = Y_mis;
  
}

model {
  
  latent_Y[1] ~ normal(stat_mu, sqrt(stat_var));
  
  latent_Y[2:N] ~ normal(lambda*mu - lambda * to_vector(latent_Y[1:(N-1)]), sigma);
  
  
  Y_obs[N_obs] ~ normal(latent_Y, epsilon);
  
  
  // Priors
  lambda ~ normal(0, 1);
  sigma  ~ normal(0, 1);
  mu     ~ normal(0, 1);
  epsilon ~ normal(0, 1);
}

generated quantities {
  
  real phi = -lambda;
  real c = lambda*mu;
  real stat_variance = (sigma^2)/(1 - lambda^2);
  
}



