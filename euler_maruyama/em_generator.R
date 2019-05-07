
drift <- function(x) {
  -0.5*(x - 0)
 
  # (4*x - 1)/x^2
  
  # 4*x*(-3 - 3*x + 2*x^2)
  
  # 0.1 - x + 1*(x^2)/(1 + x^2)
  
  # x*(1 - .5*x)*(1 + x)

  # x*(1 - .5*x)*(1 + .5*x)
}

dispersion <- function(x) {
  
  .1
  # 
  # sqrt(10*x*(1-x))
  
  # .25*abs(sinh(x)) + .25

  
  # sqrt(.10*x)
}


D_dispersion <- function(x) {
  
  1*x/sqrt(1 + x^2)

}



grid <- seq(from = 0, to = 100, by = 0.02)
# grid1 <- seq(from = 0, to = 500, by = 0.01)
# grid2 <- seq(from = 0, to = 100, by = 0.1)
# grid3 <- seq(from = 0, to = 100, by = 0.9)
# grid_list <- list(grid1 = grid1, grid2 = grid2, grid3 = grid3)
seed = sample(1:1000, 1)



## Euler - Maruyama ************************ ####
em_generator <- function(y0, grid, seed = NULL, milstein = FALSE) {
  
  if(is.null(seed)) {
    set.seed(sample(1:10000, 1))
  } else {
    set.seed(seed)
  }
        
  
      y <- y0
      for(i in 2:length(grid)) {
        
        y_prev <- y[i-1]
        delta_t <- grid[i] - grid[i-1]
        
        a <- y_prev + drift(y_prev)*delta_t + dispersion(y_prev)*rnorm(1, 0, sqrt(delta_t))
        
        if(milstein) {
          a <- a + .5*dispersion(y_prev)*D_dispersion(y_prev)*(rnorm(1, 0, delta_t) - delta_t)
        }
        
        y <- c(y, a)
        
        y
        
      }
  
      return(y)
      
}


inits <- c(1, .5)
# inits <- c(2, -2)
# # 
df <- lapply(inits, function(y0) em_generator(y0 = y0, grid = grid, seed = NULL)) %>%
  do.call(cbind, .) %>%
  as.data.frame() %>%
  cbind(time = grid) %>%
  set_colnames(c(paste0("y", 1:length(inits)), "x"))


df %>%
  melt(id.vars = "x") %>%
  ggplot(aes(x = x, y = value, color = variable)) +
  geom_line() +
  scale_color_tron() +
  labs(title = "Euler-Maruyama") +
  theme(legend.position = "none") 


# lapply(names(grid_list), function(grid) {
#   em_generator(y0 = 0, grid = grid_list[[grid]], seed = seed) %>% 
#     as.data.frame() %>% 
#     cbind(x = grid_list[[grid]], grid = grid) %>% 
#     set_colnames(c("y", "x", "grid"))
#     
# } ) %>%
#   do.call(rbind, .) %>%
#   ggplot(aes(x = x, y = y)) +
#   geom_line() +
#   scale_color_tron(labels = c("0.01", "0.05", "0.5"), name = "Step size") +
#   labs(title = "Euler-Maruyama with different step sizes") +
#   facet_wrap(~grid, ncol = 1)
#   


# lapply(inits, function(y0) em_generator(y0 = y0, grid = grid, seed = seed)) %>%
#   do.call(cbind, .) %>%
#   as.data.frame() %>%
#   cbind(lapply(inits, function(y0) em_generator(y0 = y0, grid = grid, seed = seed)) %>%
#           do.call(cbind, .) %>%
#           as.data.frame()) %>%
#   cbind(time = grid) %>%
#   set_colnames(c(paste0("em_", 1:length(inits)), paste0("m_", 1:length(inits)), "x")) %>%
#   melt(id.vars = "x") %>%
#   ggplot(aes(x = x, y = value, color = variable)) +
#   geom_line() +
#   scale_color_simpsons() +
#   labs(title = "Compare Euler-Maruyama and Milstein")

