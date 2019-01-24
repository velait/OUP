// Student t OUP

functions {
  
  // compute shape matrix 
  matrix shape(real student_df, real sigma, real lambda, real[] time, int T){
    matrix[T, T] covariance;
    matrix[T, T] shape;
    
    for(i in 1:T) {
      for(j in 1:T) {
        covariance[i, j] = ((sigma^2)./(2*lambda))*exp(-lambda*fabs(time[i]-time[j]));
      }
    }
    
    // shape = ((student_df-2)/student_df)*covariance;
    shape = covariance;
    
    return shape;
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
  int<lower=0> T;           //number of time points
  int<lower=0> N;           //number of series
  matrix[N, T] Y;           //observations
  real time[T];             //observation times
  real<lower=2> student_df;   // degrees of freedom
}

parameters {
  
  // vector [N] lambda_raw;
  // vector<lower=0, upper=1> [N] lambda;
  vector<lower=0> [N] inv_lambda;
  // vector<lower=0> [N] sigma_raw;
  vector<lower=0> [N] sigma;
  vector [N] mu_raw;
  // vector [N] mu;
  
  // real<lower=0> student_df;
  
  // hyperparameters 
  real<lower=2> inv_lambda_mean;
  // real<lower=2> lambda_mean;
  real<lower=1> inv_lambda_sd;
  real mu_mean;
  real<lower=0> mu_sd;
  real<lower=0> sigma_mean;
  real<lower=0> sigma_sd;
}

transformed parameters {
  // non-center mu and sigma
  vector [N] mu = mu_mean + mu_sd*mu_raw;
  
  vector<lower = 0> [N] lambda = rep_vector(1, N)./inv_lambda;


  // mode and variance of gamma prior
  // real<lower=0> inv_gamma_mode = lambda_sd/(lambda_mean + 1);
  // real<lower=0> inv_gamma_variance = (lambda_sd^2)/((lambda_mean-1)^2*(lambda_mean-2));
}

model {
  
  // for(i in 1:N) {
  //   
  //   Y[i, 1] ~ normal(mu[i], sqrt(2*lambda[i]*sigma[i]));
  //   
  //   
  //   for(j in 2:T) {
  //     
  //     real delta_t = time[j] - time[j-1];
  //     
  //     Y[i, j] ~ normal(mu[i] - (mu[i] - Y[i, j-1])*exp(-lambda[i]*delta_t), sqrt(2*lambda[i]*sigma[i])*(1-exp(-2*lambda[i]*delta_t)));
  //     
  //   }
    
    
    // OR
    
    for(i in 1:N) {
      Y[i] ~ multi_student_t(student_df,
      rep_vector(mu[i], T),
      ((student_df-2)/student_df)*shape(student_df, sigma[i], lambda[i], time, T));

    }
    
    
    
    // OR

  //     for(i in 1:N) {
  // 
  //   real kappa = lambda[i]*sigma[i];
  //   real log_kappa = log(kappa);
  // 
  //   target += -0.5*log_kappa - (Y[i, 1] - mu[i])/kappa;
  // 
  // 
  //   for(j in 2:T) {
  // 
  //     real delta_t = time[j] - time[j-1];
  //     real helper = 1 - exp(-2*lambda[i]*delta_t);
  // 
  // 
  //     target += -0.5*kappa - 0.5*log(helper) - (Y[i, j] - mu[i] - exp(-lambda[i]*delta_t)*(Y[i, j-1] - mu[i]))*pow(kappa*helper, -1);
  // 
  //   }
  // 
  // 
  // }
  
  // priors
  // lambda ~ inv_gamma(lambda_mean, lambda_sd);
  // lambda ~ inv_gamma(5, lambda_sd);
  inv_lambda ~ inv_gamma(inv_lambda_mean, inv_lambda_sd);
  // lambda_raw ~ normal(0, 1);
  // lambda ~ lognormal(lambda_mean, lambda_sd);
  
  
  // mu     ~ normal(mu_mean, mu_sd);
  mu_raw ~ normal(0, 1);
  
  sigma  ~ normal(sigma_mean, sigma_sd);
  // sigma_raw ~ normal(0, 1);
  // sigma ~ lognormal(sigma_mean, sigma_sd);
  
  // student_df ~ normal(0, 100);
  
  //hyper priors
  inv_lambda_mean ~ normal(5, 20);
  inv_lambda_sd ~ normal(5, 20);
  // lambda_mean ~ normal(0, 5);
  // lambda_sd ~ normal(0, 5);
  // 
  
  mu_mean ~ normal(5, 5);
  mu_sd ~ normal(0, 5);
  
  sigma_mean ~ normal(0, 5);
  sigma_sd ~ normal(0, 5);
  
  
}