library(Matrix)
library(matrixcalc)

## Parameters ************ ####
mu <- c(0, 0)

# Mean reversion
B <- matrix(c(.25, 0,
              0, 3), ncol=2, nrow=2, byrow = T)
is.positive.definite(B)

# noise
G <- matrix(c(2, 0,
              0, 1), ncol=2, nrow=2, byrow = T)
is.positive.definite(G)

is.positive.definite(B*G + G*t(B))


## Data ****************** ####

oup_2d_length <- 100

Y <- oup_2d(length = oup_2d_length, Y0 = c(0, 0), mu = mu, B = B, G = G)


Y %>% 
  cbind(time = 1:nrow(.)) %>% 
  as.data.frame() %>% 
  ggplot(aes(x = x, y = y, color = time)) + 
  geom_point(size = 1) +
  geom_path() +
  geom_point(aes(x = mu[1], y = mu[2]), color = "black", size = 2) +
  coord_cartesian(xlim = -5:5, ylim = -5:5) +
  scale_color_gradient(low = "blue", high = "red")


Y %>%
  melt() %>%
  ggplot(aes(x = Var1, y = value, color = Var2)) +
  geom_line()


# ## Stan ****************** ####
fit_2d_oup <- stan_model("2d_oup/2d_fitter.stan")


samples <- sampling(fit_2d_oup,
                    list(N = 100, x = 1:oup_2d_length, y = Y[, 1:2]),
                    iter = 1000,
                    chains = 1)


summary(samples)$summary[, c("25%", "75%")] %>%
  as.data.frame() %>%
  set_colnames(c("lower", "upper")) %>%
  rownames_to_column("parameter") %>%
  filter(parameter != "lp__") %>%
  filter(!grepl("_G", parameter)) %>%
  mutate(true_value = c(mu[1], mu[2],
                        B[1, 1], B[1, 2], B[2, 1], B[2, 2],
                        G[1, 1], G[1, 2], G[2, 1], G[2, 2])) %>%
  ggplot() +
  geom_point(aes(x = parameter, y = true_value), color = "red") +
  geom_errorbar(aes(x = parameter, ymin =lower, ymax =upper), width = 0.1)
