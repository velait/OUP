functions {
  matrix cov_exp_abs(real[] time, real alpha, real rho, int N){
    matrix[N, N] covariance;
    
    
    for(i in 1:N) {
      for(j in 1:N) {
        covariance[i, j] = square(alpha)*exp(-fabs(time[i]-time[j])/square(rho));
      }
    }
    
    return covariance;
  }
}

data {
  int<lower=1> N;
  real x[N];
  vector[N] y;
}
transformed data {
  vector[N] mu = rep_vector(0, N);
  
  real lambda1 = -log(0.01)*sqrt(2*10);
  real lambda2 = -log(0.01)/10.0;
  
}
parameters {
  real<lower=0> rho;
  real<lower=0> alpha;
  real<lower=0> sigma;
}
model {
  matrix[N, N] L_K;
  matrix[N, N] K = cov_exp_abs(x, alpha, rho, N);
  real sq_sigma = square(sigma);
  
  // diagonal elements
  for (n in 1:N) {
    K[n, n] = K[n, n] + sq_sigma;
  }
  
  
  L_K = cholesky_decompose(K);

  y ~ multi_normal_cholesky(mu, L_K);
  
  // priors
  sigma ~ normal(0, 1);
  target += log(lambda1) + log(lambda2) -(3.0/2.0)*log(rho) - lambda1*sqrt(2*rho) - lambda2*alpha;
  
  
  
}

generated quantities {
  
  real oup_sigma = sqrt(2)*alpha/rho;
  
}