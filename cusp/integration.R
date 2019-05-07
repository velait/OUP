## indefinite integration with Stan via ODE integrator and
## transformation

library(rstan)

expose_stan_functions("cusp/integrate_ode.stan")

curve(integrand(x)/(integrate(integrand, -Inf, Inf)$value), -10, 10)

theta <- c(1, 2, 1, 1)

## u = inv_logit(r)
## du = dr * d inv_logit(r) / dr

## 1/(1 + exp(-x))

## - 1/(1 + exp(-x))^2 * (-exp(-x))
## 1/(1 + exp(-x))  * ( exp(-x) )/(1+exp(-x))
## = inv_logit(x) * inv_logit(-x)

logit <- binomial()$linkfun
inv_logit <- binomial()$linkinv

trans_integrand <- function(r) {
    x <- logit(r)
    exp((theta[1]*x + 0.5*theta[2]*x^2 - 0.25*x^4)/theta[4]) * 1/(r * (1-r))
}

curve(trans_integrand, 0, 1)

## compare R naive
integrate(integrand, -Inf, Inf)
## with R transformed
integrate(trans_integrand, 0,1)
## and the Stan ode version
ode_integrate(theta)


