// data {
//   int<lower=1> N;
//   real x[N];
//   real sigma;
//   real lambda;
//   real error;
// }
// 
// transformed data {
//   real kappa = (square(sigma)/(2*lambda));
//   matrix[N, N] K;
//   vector[N] mu = rep_vector(0, N);
//   for (i in 1:(N - 1)) {
//     K[i, i] = kappa + square(error);
//     for (j in (i + 1):N) {
//       K[i, j] = kappa*exp(-0.5*lambda*fabs(x[i] - x[j]));
//       K[j, i] = K[i, j];
//     }
//   }
//   K[N, N] = kappa + square(sigma);
// }
// 
// parameters {
//   vector[N] y;
// }
// 
// model {
//   y ~ multi_normal(mu, K);
// }



data {
  int<lower=1> N;
  real x[N];

  real<lower=0> rho;
  real<lower=0> alpha;
  real<lower=0> sigma;
}

transformed data {
  matrix[N, N] cov =   cov_exp_quad(x, alpha, rho)
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


