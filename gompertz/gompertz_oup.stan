// Gompertz model; test if same problems as with cusp model persist

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
  
  vector [N_OTUs] mu;
  vector<lower = 0> [N_OTUs] lambda;
  vector<lower=0> [N_OTUs] sigma;
  
  matrix[N_OTUs, N_obs] latent_oup;   // Latent oup process
  

  // real<lower = 0> error;
}

transformed parameters {
  // real theta[N_OTUs - 1, 4];
  // theta[, 4] = epsilon;
  // for(i in 1:3) {
    //   theta[, i] = pre_theta[, i];
    // }
  // 
  // real theta[N_OTUs, 3];
  // 
  // theta[, 1] = mu;
  // theta[, 2] = lambda;
  // theta[, 3] = sigma;
  
}


model {
  // Cusp model with Shoji-Ozaki discretization ***************
    
    
    
    // for(otu in 1:(N_OTUs)) {
    //   for(t in 1:N_obs) {
    //     
    //     if(t == 1) {
    //       y[otu, t] ~ poisson_log(lambda[otu]); 
    //     } else {
    //       real prev = y[otu, t-1];
    //       // real th[3] = theta[otu, ];
    //       real delta_t = dt[t-1];
    //       // print(th);
    //       y[otu, t] ~ poisson_log(mu[otu] - (mu[otu] - prev)*exp(-lambda[otu]*delta_t),
    //                                   ((sigma[otu]^2)/(2*lambda[otu]))*(1 - exp(-2*lambda[otu]*delta_t)));
    //       
    //       // y[otu, t] ~ lognormal(latent_oup[otu, t], error);
    //     }
    //     
    //     
    //   }
    // }
  
  
  
    for(otu in 1:(N_OTUs)) {
      latent_oup[otu, 1] ~ normal(mu[otu], sigma[otu]);
      
      for(t in 2:N_obs) { 
        
        real prev = y[otu, t-1];
        real delta_t = dt[t-1];
        
        latent_oup[otu, t] ~ normal(mu[otu] - (mu[otu] - prev)*exp(-lambda[otu]*delta_t), ((sigma[otu]^2)/(2*lambda[otu]))*(1 - exp(-2*lambda[otu]*delta_t)));
        
        
        y[otu, t] ~ poisson_log(latent_oup[otu, t]);
        
        }
    }
        

  
  
  mu ~ normal(0, 2);
  lambda ~ gamma(2, 2);
  sigma ~ normal(0, 2);
  // error ~ normal(0, 2);
  
}