functions {
  
  // Exponential covariance
  matrix cov_exp(real[] time, real stat_var, real length_scale,  int T){
    
    matrix[T, T] covariance;
    
    for(i in 1:T) {
      for(j in 1:T) {
        covariance[i, j] = stat_var*exp(-fabs(time[i]-time[j])/length_scale);
      }
    }
    
    return covariance;
  }
  
}

data {
  int<lower=1> T;                 // number of timepoints
  real time[T];                   // observation times
  real<lower=0> length_scale;     // Length scale
  real<lower=0> stat_var;         // stationary variance
  real<lower=0> error;            // observation error
  int kernel;
}
transformed data {
  vector[T] zeros;
  zeros = rep_vector(0, T);
}

model {}

generated quantities {
  // real x[N];
  vector[T] y;
  vector[T] f;
  // for (n in 1:N)
  //   x[n] = uniform_rng(0,1000);
  {
    matrix[T, T] cov;
    matrix[T, T] L_cov;
    
    if(kernel == 0) {
      cov = cov_exp(time, stat_var, length_scale, T);
    } else if(kernel == 1) {
      cov = cov_exp_quad(time, sqrt(stat_var), length_scale);
    }
    
    for (t in 1:T)
      cov[t, t] = cov[t, t] + 1e-12;
    // L_cov = cholesky_decompose(cov);
    // f = multi_normal_cholesky_rng(zeros, L_cov);
    // f = multi_student_t_rng(5, zeros, cov);
  }
  for (t in 1:T)
    y[t] = normal_rng(f[t], error);
}