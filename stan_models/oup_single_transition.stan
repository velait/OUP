data {
int<lower=0> T;             //number of time points
real time[T];               //observation times
vector[T] Y;                //observations
}

transformed data {
  
  real mu = 0;
  
}

parameters {

real<lower=0> length_scale;
real<lower=0> stat_var;          
// real mu;                      // long-term mean


// vector[T] latent_process; 
// real<lower=0> epsilon;          // error standard deviation

}

transformed parameters{
  
}

model {

Y[1] ~ normal(mu, sqrt(stat_var));

for(i in 2:T) {
  
  real dt = time[i] - time[i-1];
  
  real transition_mean = mu - (mu - Y[i-1])*exp(-pow(length_scale, -1)*dt);
  real transition_var = stat_var*(1 - exp(-2*pow(length_scale, -1))*dt);
  
  Y[i] ~ normal(transition_mean, sqrt(transition_var));
  
  // Y[i] ~ normal(latent_process[i], epsilon);
  
}

// Priors
length_scale ~ gamma(2, 2);
stat_var  ~ normal(0, 1);
mu     ~ normal(0, 1);

// epsilon ~ normal(0, 1);

}



generated quantities{
  real<lower = 0> sigma = sqrt(2*stat_var/length_scale);
  real lambda = pow(length_scale, -1);
  real unit_variance = stat_var*(1 - exp(-2/length_scale));
}