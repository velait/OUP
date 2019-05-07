# Invariant density explorations
library(KScorrect)
library(reshape2)
# Given an integrated drift function, compute invariant density

potential <- function(x) x^2*(x^2 - 2)

invariant_d <- function(x, potential, var) {
  
  y <- exp(-2*potential(x)/var^2)
  
  y/sum(y)

}


x <- -200:200/100

plot(x = x, y = , type = "l")



cbind(mixture = dmixnorm(x, mean = c(-1, 1), sd = c(.1, .1), pro = c(.5, .5)), 
      model = 65*invariant_d(x, potential, .25), 
      x = x) %>% 
  as.data.frame() %>% 
  melt(id.vars = "x") %>% 
  ggplot(aes(x = x, y = value, color = variable)) + 
  geom_line()



## Potential

potential <- function(x, invariant_d) {
  
  -log(invariant_d(x))
  
}


###


plot_f <- function(a, b, func, normalize = FALSE) {
  x <- seq(a, b, length.out = 500)
  y <- eval(parse(text = func))
  
  if(normalize) {
    y <- y/sum(y)
  }
  
  p <- data.frame(x, y) %>% 
    ggplot(aes(x = x, y = y)) + 
    geom_line()
  
  p
}


# func <- ".5*( ((1-x)/(1 + exp(-(x+1)^2 + (x-1)^2))) + ((1+x)/(1 + exp((x+1)^2 - (x-1)^2))) )"
# func <- "x^2/(4*x - 1)"

pot <- "4*x^2*(x + 1)*(x - 1)"

# func <- "4*(3 + x - 6*x^2 + 2*x^3)"

# func <- paste0("exp(-", pot," )")


plot_f(-1.5, 1.5, pot)
