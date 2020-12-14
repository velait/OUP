functions {
  
  
  real drift(real x, real[] theta) {
    
    return(theta[1]*(theta[2] - x));
    
  }
  
  
  real dispersion(real x, real[] theta, real[] psi) {
    
    // Constant
    real D = theta[3];
    
    // Normal
    D = D + psi[1]/(psi[3]*sqrt(2*pi()))*exp(-(x - psi[2])^2/(2*psi[3]^2));
    
    
    return(D);
  }
  
  
  // real lamperti(real x, real[] theta) {
  //   real z = 1/(theta[4]*(1-theta[3]))*pow(x, 1 - theta[3]);
  //   return(z);
  // }  
  
  // // Below Lamperti transformed drift, L and M for the Shoji-Ozaki scheme
  // real lamperti_drift(real x, real[] theta) {
  //   return(((theta[1]*theta[2])/theta[3])*x^(-theta[4]) - (theta[1]/theta[3])*x^(1 - theta[4]) - .5*theta[3]*theta[4]*x^(theta[4] - 1));
  // }
  
  
  
}

data {
  int<lower = 1> N_obs;    // number of time points
  vector[N_obs] x;         // Observation times
  real y[N_obs];    // Observations
}

transformed data {
  vector[N_obs-1] dt = x[2:N_obs] - x[1:(N_obs-1)];
}

parameters {
  
  // Cusp
  real alpha;
  real beta;
  real lambda;
  real<lower=0> epsilon;
  real<lower=0> r;
  
  
  // Normal
  real<lower = 0> c1;
  real mu;
  real<lower = 0> sigma; // standard deviation
  
  // // Logistic
  // real c2;
  // real k;
  // real l;
  // 
}

transformed parameters {
  
  real theta[5];
  real psi[3];
  // real phi[3];
  
  theta[1] = alpha;
  theta[2] = beta;
  theta[3] = lambda;
  theta[4] = epsilon;
  theta[5] = r;
  
  
  psi[1] = c1;
  psi[2] = mu;
  psi[3] = sigma;
  
  // 
  // real phi[3];
  // phi[1] = c2;
  // phi[2] = k;
  // phi[3] = l;
}


model {
  
  y[1] ~ normal(beta, 1);
  
  
  for(t in 2:N_obs) {
    
    real y_prev = y[t-1];
    real delta_t = dt[t-1];
    
    y[t] ~ normal(y_prev + drift(y_prev, theta)*delta_t, sqrt(delta_t*dispersion(y_prev, theta, psi)));
    
  }
  
  
  alpha ~ normal(0, 2);
  beta ~ normal(0, 2);
  lambda ~ normal(0, 2);
  epsilon ~ gamma(2, 2);
  r ~ normal(2, 2);
  
  
  
  c1 ~ normal(0, 2);
  mu ~ normal(0, 2);
  sigma ~ gamma(2, 2);
  
  // 
  // real phi[3];
  // phi[1] = c2;
  // phi[2] = k;
  // phi[3] = l;
  
}