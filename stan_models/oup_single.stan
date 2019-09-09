// Student's t OUP

functions {

// compute shape matrix 
matrix covariance(real sigma, real lambda, real[] time, int T){
  
matrix[T, T] covariance;

for(i in 1:T) {
for(j in 1:T) {
covariance[i, j] = ((sigma^2)./(2*lambda))*exp(-lambda*fabs(time[i]-time[j]));
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

real<lower=0> lambda;        // mean-reversion
real<lower=0> sigma;           // volatility
real mu;                       // long-term mean

}

transformed parameters {



}

model {

// Sample observations from a multivariate t
// Y ~ multi_normal(vectorize(mu, T),
//                   covariance(sigma, lambda, time, T));


// Priors
lambda ~ inv_gamma(10, 5);
sigma  ~ normal(0.5, 1);
mu     ~ normal(0, 1);

}



