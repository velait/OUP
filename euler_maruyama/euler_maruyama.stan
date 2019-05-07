// Euler-Maruyama discretization for the cusp model 

functions {
  // real drift(real x, real alpha, real beta, real lambda, real r) {
  //   real y = r*(alpha + beta*(x - lambda) - pow(x - lambda, 3));
  //   return(y);
  // }
  
    real drift(real x, real alpha, real beta, real lambda) {
    real y = alpha + beta*(x - lambda) - pow(x - lambda, 3);
    return(y);
  }
  
}

data {
  int<lower = 1> N;
  real x[N];  
  real y[N];
}

parameters {
  
  // ordered[3]  mu;
  // vector<lower = 0> [2] mu;
  real alpha;
  real beta;
  real lambda;
  // real<lower=0> r;
  real<lower=0> epsilon;
  
  
}

model {

  // for(i in 2:N) {
  //   real delta_t = x[i] - x[i-1];
  //   y[i] ~ normal(y[i-1] + drift(y[i-1], alpha, beta, lambda, r)*delta_t, delta_t*epsilon);
  // }
  
    for(i in 2:N) {
    real delta_t = x[i] - x[i-1];
    y[i] ~ normal(y[i-1] + drift(y[i-1], alpha, beta, lambda)*delta_t, delta_t*epsilon);
  }
  
  // alpha ~ normal(0, 2.08);
  // beta ~ normal(2, 0.67);
  // lambda ~ normal(0, 1);
  // epsilon ~ normal(1, 0.1);
  // r ~ normal(1, 0.1);
  
  alpha ~ normal(1.55, 1);
  beta ~ normal(-0.13, 1);
  lambda ~ normal(0, 1);
  epsilon ~ normal(1, 1);
  // r ~ normal(1, 1);
  
}