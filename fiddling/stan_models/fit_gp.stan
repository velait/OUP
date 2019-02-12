functions {
  // vector tail_delta(vector y, vector theta, real[] x_r, int[] x_i) {
  //   vector[2] deltas;
  //   deltas[1] = inv_gamma_cdf(theta[1], exp(y[1]), exp(y[2])) - 0.01;
  //   deltas[2] = 1 - inv_gamma_cdf(theta[2], exp(y[1]), exp(y[2])) - 0.01;
  //   return deltas;
  // }
  // 
  // real[] min_delta_x(real[] x) {
  //   real diff[size(x) - 1];
  //   
  //   for(i in 2:size(x)) {
  //     diff[i] = fabs(x[i] - x[i-1]);
  //     
  //   }
  //   
  //   return(diff);
  // }
}

data {
  int<lower=1> N;
  real x[N];
  vector[N] y;
}

transformed data {
  vector[N] mu = rep_vector(0, N); // mean vector
  
  // vector[2] y_guess = [log(10), log(20)]'; // initial guess for algebra_solver
  // 
  // vector [N-1] diff =  to_vector(x[2:N]) - to_vector(x[1:(N-1)]);
  // 
  // 
  // vector[2] theta = [min(diff), max(diff)]';
  // 
  // vector[2] guess;
  // real x_r[0];
  // int x_i[0];
  // 
  // guess = algebra_solver(tail_delta, y_guess, theta, x_r, x_i);
  
}
parameters {
  real<lower=0> rho;
  real<lower=0> alpha;
  real<lower=0> sigma;
}
model {
  matrix[N, N] L_K;
  matrix[N, N] K = cov_exp_quad(x, alpha, rho);
  real sq_sigma = square(sigma);
  
  // diagonal elements
  for (n in 1:N)
    K[n, n] = K[n, n] + sq_sigma;
  
  L_K = cholesky_decompose(K);
  
  rho ~ inv_gamma(4, 10);
  alpha ~ normal(0, 2);
  sigma ~ normal(0, 1);
  
  y ~ multi_normal_cholesky(mu, L_K);
}