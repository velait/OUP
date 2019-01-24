functions {
  // compute shape matrix 
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
  
  real<lower=0> rho;
  real<lower=0> alpha;
  real<lower=0> sigma;
}

transformed data {
  matrix[N, N] cov =   cov_exp_abs(x, alpha, rho, N)
  + diag_matrix(rep_vector(1e-10, N));
  matrix[N, N] L_cov = cholesky_decompose(cov);
}

parameters {}
model {}

generated quantities {
  vector[N] f = multi_normal_cholesky_rng(rep_vector(0, N), L_cov);
  vector[N] y;
  for (n in 1:N)
    y[n] = normal_rng(f[n], sigma);
}