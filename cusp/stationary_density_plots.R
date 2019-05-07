

x_grid <- seq(from = -5, to = 5, length.out = 100)

cusp <- data.frame(x = x_grid, 
                   Unimodal = cc_density(x = x_grid, r = .5, alpha = 0, beta = -5, lambda = 0, epsilon = 2),
                   Bimodal = cc_density(x = x_grid, r = 1, alpha = -0.25, beta = 2, lambda = 0, epsilon = 1.5),
                   Skewed = cc_density(x = x_grid, r = .5, alpha = -3, beta = 2, lambda = 0, epsilon = 2),
                   Fat_tailed = cc_density(x = x_grid, r = .25, alpha = 1.5, beta = 1, lambda = 0, epsilon = 2))

ggplot(cusp %>% melt(id.vars = "x"), aes(x = x, y = value)) + 
  geom_line() + 
  labs(x = "", y = "", title = "Cusp model stationary densities", subtitle = "") +
  facet_wrap(~variable)






cusp2 <- data.frame(x = x_grid, 
                   Unimodal = cc_density(x = x_grid, r = .5, alpha = 0, beta = -5, lambda = 0, epsilon = 2),
                   Bimodal = cc_density(x = x_grid, r = 1, alpha = -0.25, beta = 2, lambda = 0, epsilon = 1.5),
                   Right_skewed = cc_density(x = x_grid, r = .4, alpha = -3, beta = 2, lambda = 0, epsilon = 2),
                   Left_skewed = cc_density(x = x_grid, r = .4, alpha = 3, beta = 2, lambda = 0, epsilon = 2))



ggplot(cusp2 %>% melt(id.vars = "x"), aes(x = x, y = value)) + 
  geom_line() + 
  labs(x = "", y = "", title = "Cusp model stationary densities", subtitle = "") +
  facet_wrap(~variable)