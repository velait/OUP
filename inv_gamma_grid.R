

x <- 1:1000/100
shape <- 
rate <- 10*log(shape)


inv_gamma_grid <- lapply(1:length(shape), function(i) {
  
  y <- dinvgamma(x, shape[i], rate[i])
  
  data.frame(y = y, x = x, i = i)
}) %>% do.call(rbind, .)


inv_gamma_grid %>% 
  ggplot(aes(x = x, y = y, color = as.factor(i))) +
  geom_line()



x <- 1:1000/100
shape <- 15
rate <- 10 + (1:10)/10


inv_gamma_grid <- lapply(1:length(rate), function(i) {
  
  y <- dinvgamma(x, shape, rate[i])
  
  data.frame(y = y, x = x, i = i)
}) %>% do.call(rbind, .)


inv_gamma_grid %>% 
  ggplot(aes(x = x, y = y, color = as.factor(i))) +
  geom_line()
