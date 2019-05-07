// Shoji approximation for the cusp model 

functions {
  
  // Drift
  real drift(real x, real[] theta) {
    real y = theta[1] + theta[2]*(x - theta[3]) - pow(x - theta[3], 3);
    return(y);
  }
  
  
  real L(real x, real[] theta) {
    return(theta[2] - 3*(x - theta[3])^2);
  }
  
  real M(real x, real[] theta) {
    return((theta[4]/2)*(- 6*(x - theta[3])));
  }
  
  real Ax(real x, real[] theta, real dt) {
    
    real l = L(x, theta);
    real m = M(x, theta);
    real d = drift(x, theta);
    real k = exp(l*dt) - 1;
    
    return(x + d*k/l + m*(k - l*dt)/l^2);
  }
  
  
  real B(real x, real[] theta, real dt) {
    
    real l = L(x, theta);
    return(sqrt(theta[4])*sqrt((exp(2*l*dt) - 1)/(2*l)));
  }
  
  
}

data {
  int<lower = 1> N;
  real x[N];  
  real y[N];
}

parameters {
  
  real theta[4];
  
  // theta[1] = alpha
  // theta[2] = beta
  // theta[3] = lambda
  // theta[4] = epsilon
}

model {
  
  y[1] ~ normal(theta[4], 1);
  
  for(i in 2:N) {
    real dt = x[i] - x[i-1];
    y[i] ~ normal(Ax(y[i-1], theta, dt), B(y[i-1], theta, dt));
  }
  
  theta[1] ~ normal(0, 2);
  theta[2] ~ normal(0, 2);
  theta[3] ~ normal(0, 2);
  theta[4] ~ normal(0, 2);
}