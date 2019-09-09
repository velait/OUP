functions {
  
  // compute shape matrix 
  matrix cov_exp(real kappa, real inv_lambda, real[] time, int T){
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
  // int<lower=1> N;
  // int<lower=1> N_pred;
  // vector[N] y;
  // real x[N];
  // real x_pred[N_pred];
  // 
  
    
  int<lower=0> T;             //number of time points
  real time[T];               //observation times
  vector[T] Y;                //observations


}
parameters {
  
  real<lower=0> inv_lambda;
  real<lower=0> kappa;
  real<lower=0> epsilon;
  vector[T] eta;
}
transformed parameters {
  vector[T] f;
  {
    matrix[T, T] L;
    matrix[T, T] K;
    K = cov_exp(kappa, inv_lambda, time, T);
    for (n in 1:T)
      K[n, n] = K[n, n] + 1e-12;
    L = cholesky_decompose(K);
    f = L * eta;
  } }
model {
  inv_lambda ~ gamma(2, .2);
  kappa ~ normal(0, 1);
  epsilon ~ normal(0, 1);
  eta ~ normal(0, 1);
  Y ~ normal(f, epsilon);
}
generated quantities {
  // vector[N_pred] f_pred;
  // vector[N_pred] y_pred;
  // f_pred = gp_pred_rng(x_pred, y, x, alpha, length_scale, sigma);
  // for (n in 1:N_pred)
    //   y_pred[n] = normal_rng(f_pred[n], sigma);
  
  real<lower = 0> sigma = sqrt(2*kappa/inv_lambda);  
  real lambda = pow(inv_lambda, -1);
  real unit_variance = kappa*(1 - exp(-2*inv_lambda));
  // vector mu = eta;
}




