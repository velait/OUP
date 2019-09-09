// Shoji-Ozaki approximation for the cusp model plus observation model

functions {
  
  // Drift
  real drift(real x, vector theta) {
    real y = theta[5]*(theta[1] + theta[2]*(x - theta[3]) - pow(x - theta[3], 3));
    return(y);
  }
  
  // Shoji-Ozaki Discretization
  real L(real x, vector theta) {
    return(theta[5]*(theta[2] - 3*(x - theta[3])^2));
  }
  
  real M(real x, vector theta) {
    return(theta[5]*((theta[4]/2)*(- 6*(x - theta[3]))));
  }
  
  real A(real x, vector theta, real dt) {
    
    real l = L(x, theta);
    real m = M(x, theta);
    real d = drift(x, theta);
    real k = exp(l*dt) - 1;
    
    return(x + d*k/l + m*(k - l*dt)/l^2);
  }
  
  real B_square(real x, vector theta, real dt) {
    
    real l = L(x, theta);
    
    // return(sqrt(theta[4])*sqrt((exp(2*l*dt) - 1)/(2*l)));
    return(theta[4]*(exp(2*l*dt) - 1)/(2*l));
  }
  

  
}

data {
  int<lower = 1> N_obs;    // number of time points
  int<lower = 1> N_OTUs;   // number of OTUs
  vector[N_obs] x;         // Observation times
  real y[N_OTUs, N_obs];    // OTU count matrix
}

transformed data {
  vector[N_obs-1] dt = x[2:N_obs] - x[1:(N_obs-1)];
}

parameters {
  
  // theta[1] = alpha
  // theta[2] = beta
  // theta[3] = lambda
  // theta[4] = epsilon
  // theta[5] = r
  
  // vector[N_OTUs] alpha_std;
  vector[N_OTUs] alpha;
  
  // vector[N_OTUs] beta_std;
  vector[N_OTUs] beta;
  
  vector [N_OTUs] lambda;
  vector<lower=0>[N_OTUs] epsilon;
  // vector<lower=0>[N_OTUs] epsilon_std;
  vector<lower=0>[N_OTUs] r;
  // vector<lower=0>[N_OTUs] r_std;
  
  // Hyperparameters
  
  real alpha_1;
  real<lower=0> alpha_2;
  real beta_1;
  real<lower=0> beta_2;
  real lambda_1;
  real<lower=0> lambda_2;
  real<lower=0> epsilon_1;
  real<lower=0> epsilon_2;
  real<lower=0> r_1;
  real<lower=0> r_2;
  
}

transformed parameters {

  // vector[N_OTUs] alpha = alpha_1 + alpha_2*alpha_std;
  // vector[N_OTUs] beta = beta_1 + beta_2*beta_std;
  // vector[N_OTUs] epsilon = epsilon_1 + epsilon_2*epsilon_std;
  // vector[N_OTUs] r = r_1 + r_2*r_std;


  matrix[N_OTUs, 5] theta;
  
  theta[, 1] = alpha;
  theta[, 2] = beta;
  theta[, 3] = lambda;
  theta[, 4] = epsilon;
  theta[, 5] = r;
  
  
  
  
  
}

model {
  
  
  for(otu in 1:(N_OTUs)) {
    
    y[otu, 1] ~ normal(lambda[otu], epsilon[otu]);
    
    for(t in 2:N_obs) {
      real prev = y[otu, t-1];
      vector[5] th = to_vector(theta[otu, ]);
      real delta_t = dt[t-1];
      
      y[otu, t] ~ normal(A(prev, th, delta_t),
                                   sqrt(B_square(prev, th, delta_t)));
    }
  }
  

  
  // alpha ~ normal(0, 1);
  alpha ~ normal(alpha_1, alpha_2);
  beta ~ normal(beta_1, beta_2);
  // beta ~ normal(0, 1);
  lambda ~ normal(lambda_1, lambda_2);
  // epsilon ~ normal(0, 1);
  epsilon ~ normal(epsilon_1, epsilon_2);
  // epsilon ~ gamma(epsilon_1, epsilon_2);
  r ~ gamma(r_1, r_2);
  // r ~ normal(0, 1);
  // r ~ normal(r_1, r_2);
  
  
  // Hyperpriors
  alpha_1 ~ normal(0, 1);
  alpha_2 ~ normal(0, 1);
  
  beta_1 ~ normal(0, 1);
  beta_2 ~ normal(0, 1);
  
  lambda_1 ~ normal(0, 1);
  lambda_2 ~ normal(0, 1);
  
  epsilon_1 ~ normal(0, 1);
  epsilon_2 ~ normal(0, 1);
  
  r_1 ~ normal(0, 5);
  r_2 ~ normal(0, 5);
}