library(Matrix)
library(matrixcalc)

## Parameters ************ ####

D <- 5

mu <- rep(0, D)

# Mean reversion
B <- matrix(c(1, 0, 0,
              0, 1, 0,
              0, 0, 1), ncol=D, nrow=D, byrow = T)
B <- diag(1, D, D)
is.positive.definite(B)

# noise
G <- matrix(c(.5, 0, 0,
              0, .5, 0,
              0, 0, .5), ncol=D, nrow=D, byrow = T)
G <- diag(1, D, D)
is.positive.definite(G)

is.positive.definite(B*G + G*t(B))

Y0 <- rep(0, D)

## Data ****************** ####

oup_multi_d_length <- 25

Y <- oup_multi_d(length = oup_multi_d_length, Y0 = Y0, mu = mu, B = B, G = G) %>% 
  cbind(time = 1:oup_multi_d_length)

# ggpairs(Y[, 1:D])



# ## Stan ****************** ####
multi_d_fitter <- stan_model("2d_oup/multi_d_fitter.stan")


samples <- sampling(multi_d_fitter,
                    list(N = oup_multi_d_length,
                         x = 1:oup_multi_d_length,
                         y = Y[, 1:D]),
                    iter = 1000,
                    chains = 1, 
                    init = list(list(Omega_G = diag(1, D, D),
                                     B = diag(1, D, D),
                                     G = diag(1, D, D))))




true_values <- c(mu,
                 B %>% t %>% as.vector,
                 G %>% t %>% as.vector)

summary(samples)$summary[, c("25%", "75%")] %>%
  as.data.frame() %>%
  set_colnames(c("lower", "upper")) %>%
  rownames_to_column("parameter") %>%
  filter(parameter != "lp__") %>%
  filter(!grepl("_G", parameter)) %>%
  mutate(true_value = true_values) %>%
  ggplot() +
  geom_point(aes(x = parameter, y = true_value), color = "red") +
  geom_errorbar(aes(x = parameter, ymin =lower, ymax =upper), width = 0.1) +
  coord_flip()
