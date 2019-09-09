


functions {
  
  
  
  matrix matern_1_5_cov(real[] x, real sigma, real rho){
    
    int T = size(x);
    
    matrix[T, T] covariance;
    
    for(i in 1:T) {
      for(j in 1:T) {
        real term = (sqrt(3)*fabs(x[i]-x[j]))/rho;
        covariance[i, j] = square(sigma)*(1 + term)*exp(-term);
      }
    }
    
    return covariance;
  }
  
  matrix matern_1_5_cov2(real[] x1, real[] x2, real sigma, real rho){
    
    int T1 = size(x1);
    int T2 = size(x2);
    
    matrix[T1, T2] covariance;
    
    for(i in 1:T1) {
      for(j in 1:T2) {
       real term = (sqrt(3)*fabs(x1[i]-x2[j]))/rho;
       covariance[i, j] = square(sigma)*(1 + term)*exp(-term);
      }
    }
    
    return covariance;
  }
  
  
  
  
  
  
  
  vector gp_pred_rng(real[] x_pred,
                     vector y_is,
                     real[] x_is,
                     real sigma,
                     real rho,
                     real epsilon) {
    vector[size(x_pred)] f_pred;
    int N_pred;
    int N;
    N_pred = size(x_pred);
    N = rows(y_is);
    {
      matrix[N, N] L_Sigma;
      vector[N] K_div_y_is;
      matrix[N, N_pred] k_x_is_x_pred;
      matrix[N, N_pred] v_pred;
      vector[N_pred] f_pred_mu;
      matrix[N_pred, N_pred] cov_f_pred;
      matrix[N_pred, N_pred] nug_pred;
      matrix[N, N] Sigma;
      Sigma = matern_1_5_cov(x_is, sigma, rho);
      for (n in 1:N)
        Sigma[n, n] = Sigma[n,n] + square(epsilon);
      L_Sigma = cholesky_decompose(Sigma);
      K_div_y_is = mdivide_left_tri_low(L_Sigma, y_is);
      K_div_y_is = mdivide_right_tri_low(K_div_y_is',L_Sigma)';
                                         k_x_is_x_pred = matern_1_5_cov2(x_is, x_pred, sigma, rho);
                                         f_pred_mu = (k_x_is_x_pred' * K_div_y_is);
      v_pred = mdivide_left_tri_low(L_Sigma, k_x_is_x_pred);
      cov_f_pred = matern_1_5_cov(x_pred, sigma, rho) - v_pred' * v_pred;
                                                      nug_pred = diag_matrix(rep_vector(1e-12,N_pred));
                                                      f_pred = multi_normal_rng(f_pred_mu, cov_f_pred + nug_pred);
    }
    return f_pred;
  }
}
data {
  int<lower=1> N;
  int<lower=1> N_pred;
  vector[N] y;
  real x[N];
  real x_pred[N_pred];
}

transformed data {
  
  real<lower=0> epsilon = 0.1;
  
}

parameters {
  
  real<lower=0> rho;
  real<lower=0> sigma;
  // real<lower=0> epsilon;
  vector[N] eta;
}
transformed parameters {
  vector[N] f;
  {
    matrix[N, N] L;
    matrix[N, N] K;
    K = matern_1_5_cov(x, sigma, rho);
    for (n in 1:N)
      K[n, n] = K[n, n] + 1e-12;
    L = cholesky_decompose(K);
    f = L * eta;
  } }
model {
  rho ~ gamma(2, .5);
  sigma ~ normal(0, 1);
  // epsilon ~ normal(0, 1);
  eta ~ normal(0, 1);
  y ~ normal(f, epsilon);
}
generated quantities {
  vector[N_pred] f_pred;
  vector[N_pred] y_pred;
  f_pred = gp_pred_rng(x_pred, y, x, sigma, rho, epsilon);
  for (n in 1:N_pred)
    y_pred[n] = normal_rng(f_pred[n], epsilon);
}