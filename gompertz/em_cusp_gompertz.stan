// Euler approximation for the cusp gompertz model 

functions {
  
    // Drift
    real drift(real x, real[] theta) {
    real y = theta[1] + theta[2]*(x - theta[3]) - pow(x - theta[3], 3);
    return(y);
    }
  
  
  real dispersion(real x, real[] theta) {
    
    return(theta[4]);
    
  }
  
  
  
}

data {
  int<lower = 1> N_obs;    // number of time points
  int<lower = 1> N_OTUs;   // number of OTUs
  vector[N_obs] x;         // Observation times
  int y[N_OTUs, N_obs];    // OTU count matrix
}

transformed data {
  vector[N_obs-1] dt = x[2:N_obs] - x[1:(N_obs-1)];
}

parameters {
  
  // real theta[N_OTUs, 4];
  
  // real pre_theta[N_OTUs-1, 3];
  // real<lower = 0> epsilon[N_OTUs-1];
  // theta[1] = alpha
  // theta[2] = beta
  // theta[3] = lambda
  // theta[4] = epsilon
  
  real alpha[N_OTUs];
  real beta[N_OTUs];
  real lambda[N_OTUs];
  real<lower=0> epsilon[N_OTUs];
  real<lower = 0> r[N_OTUs];
  
  matrix[N_OTUs, N_obs] latent_cusp;   // Latent cusp process
  
}

transformed parameters {
  // real theta[N_OTUs - 1, 4];
  // theta[, 4] = epsilon;
  // for(i in 1:3) {
    //   theta[, i] = pre_theta[, i];
    // }
  
  real theta[N_OTUs, 5];
  
  theta[, 1] = alpha;
  theta[, 2] = beta;
  theta[, 3] = lambda;
  theta[, 4] = epsilon;
  theta[, 5] = r;
}


model {
  // Cusp model with Shoji-Ozaki discretization ***************
    
  
  for(otu in 1:N_OTUs) {
    
    latent_cusp[otu, 1] ~ normal(lambda[otu], 1);
    
    for(t in 2:N_obs) {
      
      real y_prev = y[otu, t-1];
      real delta_t = dt[t-1];
      real th[5] = theta[otu, ];
      
      
      latent_cusp[otu, t] ~ normal(y_prev + drift(y_prev, th)*delta_t, sqrt(delta_t*dispersion(y_prev, th)));
      
      y[otu, t] ~ poisson_log(latent_cusp[otu, t]);
       
    }
  
    
  }  
  
  alpha ~ normal(0, 2);
  beta ~ normal(0, 2);
  lambda ~ normal(0, 2);
  epsilon ~ normal(0, 2);
  r ~ normal(0, 2);
  
  
  
}