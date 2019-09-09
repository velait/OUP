// Shoji approximation for the cusp model 

functions {
  
  // Drift
  real drift(real x, real[] theta) {
    real y = theta[5]*(theta[1] + theta[2]*(x - theta[3]) - pow(x - theta[3], 3));
    return(y);
  }
  
  // Shoji-Ozaki Discretization
  real L(real x, real[] theta) {
    return(theta[5]*(theta[2] - 3*(x - theta[3])^2));
  }
  
  real M(real x, real[] theta) {
    return(-3*theta[4]*theta[5]*(x - theta[3]));
  }
  
  
  // in the model block y ~ normal(Ax, B^2), so here A = Ax
  // real A(real x, real[] theta, real dt) {
  //   
  //   real l = L(x, theta);
  //   real m = M(x, theta);
  //   real d = drift(x, theta);
  //   real k = exp(l*dt) - 1;
  //   
  //   return(x + d*k/l + m*(k - l*dt)/l^2);
  // }
  
  
  real A(real x, real[] theta, real dt) {
    
    real l = L(x, theta);
    
    real n1 = theta[1] + theta[2]*(x - theta[3]) - (x - theta[3])^3;
    real d1 = theta[2] - 3*(x - theta[3])^2;
    real n2 = 3*theta[4]*(x - theta[3])*(exp(l*dt) - 1 - l*dt);
    real d2 = theta[5]*(theta[2] - 3*(x - theta[3])^2);
    
    return(x + n1/d1 - n2/d2);
  }
  
  real B_square(real x, real[] theta, real dt) {
    
    real l = L(x, theta);
    
    // return(sqrt(theta[4])*sqrt((exp(2*l*dt) - 1)/(2*l)));
    return(theta[4]*(exp(2*l*dt) - 1)/(2*l));
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
  
  
  
  for(otu in 1:(N_OTUs)) {
    for(t in 1:N_obs) {
      
      if(t == 1) {
        y[otu, t] ~ poisson_log(latent_cusp[otu, t]); 
      } else {
        real prev = y[otu, t-1];
        real th[5] = theta[otu, ];
        real delta_t = dt[t-1];
        // print(th);
        latent_cusp[otu, t] ~ normal(A(prev, th, delta_t),
                                     sqrt(B_square(prev, th, delta_t)));
                                     
        y[otu, t] ~ poisson_log(latent_cusp[otu, t]);
      }
      

    }
  }

  alpha ~ normal(0, 2);
  beta ~ normal(0, 2);
  lambda ~ normal(0, 2);
  epsilon ~ normal(0, 2);
  r ~ normal(0, 2);
}