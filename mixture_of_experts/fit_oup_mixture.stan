
data {
  int<lower = 1> N;
  real x[N];  
  real y[N];
}

parameters {
  
  // parameters
  ordered [2] mu;
  vector<lower = 0> [2] sigma;
  vector<lower = 0> [2] lambda;
  
  // latent process; mixture weight
  vector<lower = 0, upper = 1> [N] mix_weight;
  
}

model {
  
    
  
  
target += log_mix(mix_weight[1],
              normal_lpdf(y[1] | mu[1], (square(sigma[1]) ./ (2*lambda[1]))),
              normal_lpdf(y[2] | mu[2], (square(sigma[2]) ./ (2*lambda[2]))));


  for (i in 2:N) {
    
    real delta_t = x[i] - x[i-1];
    
    target += log_mix(mix_weight[i],
                    normal_lpdf(y[i] | mu[1] - (mu[1] - y[i-1]) .* exp(-lambda[1]*delta_t), (square(sigma[1]) ./ (2*lambda[1])) .* (1 - exp(-2*lambda[1]*delta_t))),
                    normal_lpdf(y[i] | mu[2] - (mu[2] - y[i-1]) .* exp(-lambda[2]*delta_t), (square(sigma[2]) ./ (2*lambda[2])) .* (1 - exp(-2*lambda[2]*delta_t))));
                    
      // Random walk prior for the latent process
      mix_weight[i] ~ normal(mix_weight[i-1], 0.25*delta_t);
      target += log_mix(mix_weight[i-1], normal_lpdf(mix_weight[i] | 1 , .25), normal_lpdf(mix_weight[i] | 0 , .25));
      
}


// priors

mu[1] ~ normal(0, 1);
mu[2] ~ normal(1, 1);
sigma ~ normal(0, 1);
// sigma ~ beta(.5, .5);
lambda ~ inv_gamma(4, 10);
// lambda ~ normal(0, 1);
// lambda ~ normal(1, 1);





}