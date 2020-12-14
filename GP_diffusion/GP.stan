functions {
  
  // Custom function for prediction generation
  vector gp_pred_rng(real[] x_pred,
                     vector y_is,
                     real[] x_is,
                     real alpha,
                     real length_scale,
                     real sigma) {
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
      Sigma = cov_exp_quad(x_is, alpha, length_scale);
      for (n in 1:N)
        Sigma[n, n] = Sigma[n,n] + square(sigma);
      L_Sigma = cholesky_decompose(Sigma);
      K_div_y_is = mdivide_left_tri_low(L_Sigma, y_is);
      K_div_y_is = mdivide_right_tri_low(K_div_y_is',L_Sigma)';
                                         k_x_is_x_pred = cov_exp_quad(x_is, x_pred, alpha, length_scale);
                                         f_pred_mu = (k_x_is_x_pred' * K_div_y_is);
                                                      v_pred = mdivide_left_tri_low(L_Sigma, k_x_is_x_pred);
                                                      cov_f_pred = cov_exp_quad(x_pred, alpha, length_scale) - v_pred' * v_pred;
                                                      nug_pred = diag_matrix(rep_vector(1e-12,N_pred));
                                                      f_pred = multi_normal_rng(f_pred_mu, cov_f_pred + nug_pred);
    }
    return f_pred;
  }

}
data {
  int<lower=1> N;        // Number of observations
  int<lower=1> N_pred;   // Number of predicted values
  vector[N] y;           // Data: y values
  real x[N];             // Data: x values
  real x_pred[N_pred];   // x values where predictions made
  
  // Hyperparameters
  real<lower=0> length_scale;     // Length scale
  real<lower=0> stat_var;         // stationary variance
  real<lower=0> error;            // observation error

}

parameters {
  
  vector[N] eta;                  // Implicit non-centered parameterization. See Stan manual for more info
}

transformed parameters {
  
  // Generate covariance matrix
  vector[N] f;
  {
    matrix[N, N] L;
    matrix[N, N] K;
    
    K = cov_exp_quad(x, sqrt(stat_var), length_scale);
  
  // Add nugget to diagonal, for computational stability
    for (t in 1:N)
      K[t, t] = K[t, t] + 1e-12;
      
  // Cholesky transform
    L = cholesky_decompose(K);
    f = L * eta;
  }
}


model {
  
  // Latent variable generated from N(0, 1)
  eta ~ normal(0, 1);
  
  // Data is a noisy observation of latent function f
  y ~ normal(f, error);
  
}

generated quantities {
  
  // Generate predictions
  vector[N_pred] f_pred;
  vector[N_pred] y_pred;
  f_pred = gp_pred_rng(x_pred, y, x, sqrt(stat_var), length_scale, error);
  for (n in 1:N_pred)
    y_pred[n] = normal_rng(f_pred[n], error);
}