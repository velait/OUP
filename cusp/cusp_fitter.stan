functions {
  
  // real[] integrand_ode(real r, real[] f, real[] theta, real[] x_r, int[] x_i) {
  //   real df_dx[1];
  //   real x = logit(r);
  //   df_dx[1] = exp((theta[1]*(x-theta[3]) + 0.5*theta[2]*square(x-theta[3]) - 0.25*square(square(x-theta[3])))/theta[4]);
  //   df_dx[1] = 1/(r * (1-r));
  //   return(df_dx);
  // }

  real[] integrand_ode(real r, real[] f, real[] theta, real[] x_r, int[] x_i) {
    real df_dx[1];
    real x = logit(r);
    df_dx[1] = exp((theta[1]*x + 0.5*theta[2]*square(x) - 0.25*square(square(x)))/1);
    df_dx[1] = 1/(r * (1-r));
    return(df_dx);
  }

  real ode_integrate(real[] theta) {
    
    int x_i[0];
    return(
      integrate_ode_rk45(
        integrand_ode,
        rep_array(0.0, 1),
        1E-5,
        rep_array(1.0-1E-5, 1),
        theta,
        rep_array(0.0, 0),
        x_i)[1,1]
      );
      
  }
  
  
  
  // real invariant_density_lpdf(real y, real[] theta) {
  //   
  //   real log_y;
  //   log_y = (theta[1]*(y - theta[3]) + 0.5*theta[2]*pow(y - theta[3], 2) - 0.25*pow(y - theta[3], 4))/theta[4];
  //   
  //   log_y = log_y - ode_integrate(theta);
  //   
  //   return(log_y);
  // }
 
   real invariant_density_lpdf(real y, real[] theta) {
    
    real log_y;
    log_y = (theta[1]*y + 0.5*theta[2]*pow(y, 2) - 0.25*pow(y, 4))/1;
    
    log_y = log_y - log(ode_integrate(theta));
    
    return(log_y);
  }
  
  
}

data {
  int N;
  vector[N] y;
}


// transformed data {
//   real delta = 1;
// }


parameters {
  
  
  real theta[2]; // {alpha, beta, lambda, delta}
  
  // real lambda;
  // // real<lower = 0> epsilon;
  // // real<lower = 0> r;
  // real<lower = 0> delta;
  
  // real mu;
  // real<lower=0> std;
}

model {
  
  // for(i in 1:N) {
  //   target += invariant_density_lpdf(y[i] | alpha, beta, lambda, epsilon, r);
  // }
  
  //   for(i in 1:N) {
  //     target += invariant_density_lpdf(y[i] | alpha, beta, lambda, delta);
  // }
  
    for(i in 1:N) {
      target += invariant_density_lpdf(y[i] | theta);
  }
  
  // for(i in 1:N) {
  //   target += myNormal_lpdf(y[i] |mu, std);
  // }

  

  theta[1:2] ~ normal(0, 1);
  // theta[4] ~ gamma(1.05, 0.01);

  // alpha ~ normal(0, 2);
  // beta ~ normal(5, 2);
  // lambda ~ normal(5, .1);
  // delta ~ normal(5, .1);
  // epsilon ~ normal(0, 1);
  // r ~ normal(0, 2);
  
  
  // mu ~ normal(0, 1);
  // std ~ normal(0, 1);

}