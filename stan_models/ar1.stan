data {
  int<lower=0> T;             //number of time points
  real time[T];               //observation times
  vector[T] Y;                //observations
}

transformed data {
  real epsilon = 0.1; // error standard deviation
}

parameters {
  
  real mu;
  real<lower=-1, upper=1> lambda;        // mean-reversion
  real<lower=0> sigma;          
  
  
  vector[T] latent_Y;
}

transformed parameters{
  real stat_mu = (lambda*mu)/(1 + lambda);
  real<lower=0> stat_var = (sigma^2)/(1 - lambda^2);
}

model {
  
  latent_Y[1] ~ normal(stat_mu, sqrt(stat_var));

  latent_Y[2:T] ~ normal(lambda*mu - lambda*Y[1:(T-1)], sigma);
    
    
  Y ~ normal(latent_Y, epsilon);
  
  // Priors
  lambda ~ normal(0, 1);
  sigma  ~ normal(0, 1);
  mu     ~ normal(0, 1);
  
}

generated quantities {
  
  real phi = -lambda;
  real c = lambda*mu;
  real stat_variance = (sigma^2)/(1 - lambda^2);
  

}



