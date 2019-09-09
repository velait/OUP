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
  int<lower=1> T;                 // number of timepoints
  real x[T];                   // observation times
  real<lower=0> rho;              // lengthscale
  real<lower=0> alpha;         //  standard deviations
  
  real<lower=0> error;            // observation error
  
}
transformed data {
  vector[T] zeros;
  zeros = rep_vector(0, T);
}

model {}

generated quantities {
  
  vector[T] y;
  vector[T] f;
  
    {
      matrix[T, T] cov;
      matrix[T, T] L_cov;
      
      // cov = cov_exp_quad(x, alpha, rho);
      cov = cov_exp(alpha, rho, x, T);
      
      for (t in 1:T)
        cov[t, t] = cov[t, t] + 1e-12;
        
        
      L_cov = cholesky_decompose(cov);
      f = multi_normal_cholesky_rng(zeros, L_cov);
      
    }
    for (t in 1:T)
      y[t] = normal_rng(f[t], error);
}