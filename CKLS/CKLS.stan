functions {
  
  real lamperti(real x, real[] theta) {
    real z = 1/(theta[4]*(1-theta[3]))*pow(x, 1 - theta[3]);
    return(z);
  }  
  
  // Below Lamperti transformed drift, L and M for the Shoji-Ozaki scheme
  real drift(real x, real[] theta) {
    return(((theta[1]*theta[2])/theta[3])*x^(-theta[4]) - (theta[1]/theta[3])*x^(1 - theta[4]) - .5*theta[3]*theta[4]*x^(theta[4] - 1));
  }
  
  real L(real x, real[] theta) {
    return((-(theta[1]*theta[2]*theta[4])/theta[3])*x^(-theta[4] - 1) + 
    (theta[1]*(theta[4] - 1)/theta[3])*x^(-theta[4]) + 
    .5*theta[3]*theta[4]*(1-theta[4])*x^(theta[4] - 2));
  }
  
  real M(real x, real[] theta) {
    return(.5*( (theta[1]*theta[2]*theta[4]*(1+theta[4]))/theta[3]*x^(-theta[4] - 2) +
              (theta[1]*theta[4]*(1-theta[4]))/theta[3]*x^(-theta[4] - 1) +
              .5*theta[3]*theta[4]*(1-theta[4])*(theta[4]-2)*x^(theta[4] - 3)
    ));
  }
  
  
  real A(real x, real[] theta, real dt) {
    
    real l = L(x, theta);
    real m = M(x, theta);
    real d = drift(x, theta);
    real k = exp(l*dt) - 1;
    
    return(x + d*k/l + m*(k - l*dt)/l^2);
  }
  
  real B_square(real x, real[] theta, real dt) {
    
    real l = L(x, theta);
    
    return(theta[4]*(exp(2*l*dt) - 1)/(2*l));
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
  
  real z[N_obs];
  for(t in 1:N_obs) {
    
    z[t] = lamperti(y[t], theta);
    
  }
  
  
  
  
    
  for(t in 2:N_obs) {
  real prev = z[t-1];
  real delta_t = dt[t-1];
  z[t] ~ normal(A(prev, theta, delta_t),
                               sqrt(B_square(prev, theta, delta_t)));
}

  
    
  
  alpha ~ normal(0, );
  beta ~ normal(0, 21);
  epsilon ~ normal(0, 2);
  gamma ~ normal(0, 2);
  
}