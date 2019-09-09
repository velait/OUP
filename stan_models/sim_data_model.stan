functions {
  
    // compute shape matrix 
  matrix cov_exp(real[] time, real kappa, real inv_lambda,  int T){
    matrix[T, T] covariance;
    
    
    for(i in 1:T) {
      for(j in 1:T) {
        covariance[i, j] = kappa*exp(-fabs(time[i]-time[j])/inv_lambda);
      }
    }
    
    return covariance;
  }
  
}

data {
  int<lower=1> N;
  real<lower=0> length_scale;
  real<lower=0> alpha;
  real<lower=0> sigma;
  
}
transformed data {
  vector[N] zeros;
  zeros = rep_vector(0, N);
}
model {}
generated quantities {
  real x[N];
  vector[N] y;
  vector[N] f;
  for (n in 1:N)
    x[n] = uniform_rng(0,1000);
  {
    matrix[N, N] cov;
    matrix[N, N] L_cov;
    // cov = cov_exp(x, alpha, length_scale, N);
    cov = cov_exp_quad(x, alpha, length_scale);
    for (n in 1:N)
      cov[n, n] = cov[n, n] + 1e-12;
    L_cov = cholesky_decompose(cov);
    f = multi_normal_cholesky_rng(zeros, L_cov);
  }
  for (n in 1:N)
    y[n] = normal_rng(f[n], sigma);
}