

functions {
  
  // compute shape matrix 
  matrix exp_cov(real nu, real kappa, real inv_lambda, real[] time, int T){
    matrix[T, T] covariance;
    
    
    for(i in 1:T) {
      for(j in 1:T) {
        covariance[i, j] = kappa*exp(-fabs(time[i]-time[j])/inv_lambda);
      }
    }
    
    covariance = ((nu-2)/nu)*covariance;
    
    return covariance;
  }
  
  
    // compute shape matrix 
  matrix exp_cov_sq(real nu, real kappa, real inv_lambda, real[] time, int T){
    matrix[T, T] covariance;
    
    
    for(i in 1:T) {
      for(j in 1:T) {
        covariance[i, j] = square(kappa)*exp(-fabs(time[i]-time[j])/(2*inv_lambda^2));
      }
    }
    
    covariance = ((nu-2)/nu)*covariance;
    
    return covariance;
  }
  
  
  // make vector with equal inputs
  vector vectorize(real mu, int T) {
    vector[T] vec;
    for(i in 1:T) {
      vec[i] = mu;
    }
    return(vec);
  }
  
}

data {
  int<lower=0> T;             //number of time points
  real time[T];               //observation times
  vector[T] Y;                //observations
  
  // real kappa_prior[2];
  
}

transformed data {
  
  real nu = 5;
  
}

parameters {
  
  real log_inv_lambda;        // mean-reversion
  real log_kappa;           // volatility
  real mu;                       // long-term mean
  
  real<lower=0> epsilon;
  
  
  
}

transformed parameters {
  
  real inv_lambda = exp(log_inv_lambda);
  real kappa = exp(log_kappa);
  
  
  
  
}

model {
  
  // Jacobian
  target += log_inv_lambda;
  target += log_kappa;
  
  
  // Sample observations from a multivariate t
  Y ~ multi_student_t(nu,
                      rep_vector(mu, T),
                      exp_cov_sq(nu, sqrt(kappa), inv_lambda, time, T) + diag_matrix(rep_vector(square(epsilon), T)));
  
  
  // Priors
  inv_lambda ~ gamma(10, 2);
  kappa  ~ gamma(5, 20);
  mu     ~ normal(0, 1);
  
  nu ~ gamma(2, .1);
  
  epsilon ~ normal(0, 1);
  
}

generated quantities {
  real stat_var = kappa;
  real length_scale = inv_lambda;
  real lambda = pow(inv_lambda, -1);
  real sigma = sqrt(2*lambda*kappa);
  
}