# Bimodal simulator

mixture <- .4

mix_dnorm <- function(x) {
  mixture*dnorm(x, -1, .5) + (1 - mixture)*dnorm(x, 1, .5)
} 


inv_tau <- function(x) {
  
  x <- bimodal_quantile_function(x)
  
  qnorm(x)
  
}