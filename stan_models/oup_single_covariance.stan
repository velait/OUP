

functions {

// compute shape matrix 
matrix exp_cov(real kappa, real inv_lambda, real[] time, int T){
matrix[T, T] covariance;


for(i in 1:T) {
for(j in 1:T) {
covariance[i, j] = kappa*exp(-fabs(time[i]-time[j])/inv_lambda);
}
}

return covariance;
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

}

parameters {

real<lower=0> inv_lambda;        // mean-reversion
real<lower=0> kappa;           // volatility
real mu;                       // long-term mean

real<lower=0> epsilon;

}

transformed parameters {



}

model {

// Sample observations from a multivariate t
Y ~ multi_normal(rep_vector(mu, T), exp_cov(kappa, inv_lambda, time, T) + diag_matrix(rep_vector(square(epsilon), T)));


// Priors
inv_lambda ~ gamma(10, 2);
kappa  ~ gamma(5, 20);
mu     ~ normal(0, 1);

epsilon ~ normal(0, 1);

}

generated quantities {
  
  real lambda = pow(inv_lambda, -1);
  real sigma = sqrt(2*lambda*kappa);
  
}