// Student's t OUP

functions {
  
  // compute shape matrix 
  matrix shape(real student_df, real sigma, real lambda, real[] time, int T){
    matrix[T, T] covariance;
    matrix[T, T] shape;
    
    for(i in 1:T) {
      for(j in 1:T) {
        covariance[i, j] = ((sigma^2)./(2*lambda))*exp(-lambda*fabs(time[i]-time[j]));
      }
    }
    
    // shape = ((student_df-2)/student_df)*covariance;
    shape = covariance;
    
    return shape;
  }
  
  // make vector with equal inputs
  vector vectorize(real mu, int T) {
    vector[T] vec;
    for(i in 1:T) {
      vec[i] = mu;
    }
    return(vec);
  }
  
}

data {
  int<lower=0> T;             //number of time points
  real time[T];               //observation times
  vector[T] Y;                //observations
  real<lower=2> student_df;   // degrees of freedom
}

parameters {
  
  real<lower=0> lambda;        // mean-reversion
  real<lower=0> sigma;           // volatility
  real mu;                       // long-term mean
  
}

transformed parameters {
  
  
  
}

model {
  
    // Sample observations from a multivariate t
    Y ~ multi_student_t(student_df,
                        vectorize(mu, T),
                        ((student_df-2)/student_df)*shape(student_df, sigma, lambda, time, T));
    
    
    // Priors
    lambda ~ inv_gamma(10, 5);
    sigma  ~ normal(0.5, 1);
    mu     ~ normal(0, 1);
    
}
  


