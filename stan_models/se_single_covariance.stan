functions {
  
  // compute shape matrix 
  matrix cov_exp(real var_par, real lengthscale, real[] time, int T){
    matrix[T, T] covariance;
    
    
    for(j in 1:T) {
      for(i in 1:T) {
        covariance[i, j] = square(var_par)*exp(-fabs(time[i]-time[j])/square(lengthscale));
      }
    }
  
    return covariance;
  }
  
}


data {
  int<lower=1> N;           //number of observation
  real x[N];                
  vector[N] y;                //observations
}

transformed data {
  vector[N] mu = rep_vector(0, N);
}

parameters {
  
  real<lower=0> rho;        // lengthscale
  real<lower=0> alpha;      // standard deviation
  real<lower=0> sigma;      // measurement error std
  
}

transformed parameters {
  
  
  
}

model {
  
  matrix[N, N] L_K;
  // matrix[N, N] K = cov_exp_quad(x, alpha, rho);
  matrix[N, N] K =  cov_exp(alpha, rho, x, N);
  
  
  real sq_sigma = square(sigma);

  // diagonal elements
  for (n in 1:N)
    K[n, n] = K[n, n] + sq_sigma;

  L_K = cholesky_decompose(K);

  rho ~ gamma(5, 2);
  alpha ~ gamma(2, 2);
  sigma ~ normal(0, 1);

  y ~ multi_normal_cholesky(mu, L_K);
  
}

generated quantities {
  
  real length_scale = rho;
  real stat_var = square(alpha);
}
