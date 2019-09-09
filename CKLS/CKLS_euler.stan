functions {
  
  
  real drift(real x, real[] theta) {
    
    return(theta[1]*(theta[2] - x));
    
  }
  
  
  real dispersion(real x, real[] theta) {
    
    return(theta[3]*x^theta[4]);
    
  }
  
  
  real lamperti(real x, real[] theta) {
    real z = 1/(theta[4]*(1-theta[3]))*pow(x, 1 - theta[3]);
    return(z);
  }  
  
  // Below Lamperti transformed drift, L and M for the Shoji-Ozaki scheme
  real lamperti_drift(real x, real[] theta) {
    return(((theta[1]*theta[2])/theta[3])*x^(-theta[4]) - (theta[1]/theta[3])*x^(1 - theta[4]) - .5*theta[3]*theta[4]*x^(theta[4] - 1));
  }
  
  
  
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
  real alpha;
  real beta;
  real<lower=0> epsilon;
  real<lower=0> gamma;
  
}

transformed parameters {
  
  real theta[4];
  
  theta[1] = alpha;
  theta[2] = beta;
  theta[3] = epsilon;
  theta[4] = gamma;
  
}


model {
  
  y[1] ~ normal(beta, 1);
  
  
  for(t in 2:N_obs) {
    
    real y_prev = y[t-1];
    real delta_t = dt[t-1];
    
    y[t] ~ normal(y_prev + drift(y_prev, theta)*delta_t, sqrt(delta_t*dispersion(y_prev, theta)));
    
  }
  
  
  alpha ~ normal(0, 2);
  beta ~ normal(0, 2);
  epsilon ~ gamma(2, 2);
  gamma ~ normal(2, 2);
  
}