data {
  int<lower=0> T;             //number of time points
  real time[T];               //observation times
  vector[T] Y;                //observations
  
  real stat_var_approx;
  real length_scale_approx[2];
}

transformed data {
  
  real mu = 0;
  
  real length_scale_a = square(length_scale_approx[1])/(length_scale_approx[2] + .5);
  real length_scale_b = length_scale_approx[1]/(length_scale_approx[2] + .5);
  
  
}

parameters {
  
  real<lower=0> length_scale;
  real<lower=0> stat_var;          
  // real mu;                      // long-term mean
  
  real<lower=0> epsilon;
  
}

transformed parameters{
  
}

model {
  
  Y[1] ~ normal(mu, stat_var);
  
  for(i in 2:T) {
    
    real transition_mean = mu - (mu - Y[i-1])*exp(-pow(length_scale, -1));
    real transition_var = stat_var*(1 - exp(-2*pow(length_scale, -1)));
    
    Y[i] ~ normal(transition_mean, sqrt(transition_var) + epsilon);
    
  }
  
  // Priors
  length_scale ~ gamma(length_scale_a, length_scale_b);
  stat_var  ~ normal(stat_var_approx, 0.1);
  // mu     ~ normal(0, 1);
  
  epsilon ~ normal(0, 1);
  
}



generated quantities{
  real<lower = 0> sigma = sqrt(2*stat_var/length_scale);
  real lambda = pow(length_scale, -1);
  real unit_variance = stat_var*(1 - exp(-2/length_scale));
}