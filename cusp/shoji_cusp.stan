// Shoji approximation for the cusp model 

functions {
    
    // Drift
    real drift(real x, real[] theta) {
    real y = (theta[1] + theta[2]*(x - theta[3]) - pow(x - theta[3], 3));
    return(y);
    }
  
    // Shoji-Ozaki Discretization
    real L(real x, real[] theta) {
      return(theta[2] - 3*(x - theta[3])^2);
    }
  
    real M(real x, real[] theta) {
      return((theta[4]/2)*(- 6*(x - theta[3])));
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
    
    // return(sqrt(theta[4])*sqrt((exp(2*l*dt) - 1)/(2*l)));
    return(theta[4]*(exp(2*l*dt) - 1)/(2*l));
    }

    // Bijective softmax  
    // vector softmax_id(vector alpha) {
    //    vector[num_elements(alpha) + 1] alphac1;
    //    for (k in 1:num_elements(alpha))
    //      alphac1[k] = alpha[k];
    //    alphac1[num_elements(alphac1)] = 0;
    //    return softmax(alphac1);
    //   }

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
  

}

transformed parameters {
  // real theta[N_OTUs - 1, 4];
  // theta[, 4] = epsilon;
  // for(i in 1:3) {
  //   theta[, i] = pre_theta[, i];
  // }
  
  real theta[N_OTUs, 4];
  
  theta[, 1] = alpha;
  theta[, 2] = beta;
  theta[, 3] = lambda;
  theta[, 4] = epsilon;
}


model {
  // Cusp model with Shoji-Ozaki discretization ***************
  
  y[, 1] ~ normal(0, 1); 
  
  for(otu in 1:(N_OTUs)) {
    for(t in 2:N_obs) {
      real prev = y[otu, t-1];
      real th[4] = theta[otu, ];
      real delta_t = dt[t-1];
      // print(th);
      y[otu, t] ~ normal(A(prev, th, delta_t),
                                   sqrt(B_square(prev, th, delta_t)));
    }
  }

  alpha ~ normal(0, 2);
  beta ~ normal(0, 2);
  lambda ~ normal(0, 2);
  epsilon ~ normal(0, 2);
  // epsilon ~ gamma(2, 2);
}