## Simulate data from Lotka Volterra or Faust's Toy Model

## Data ******************************************** ####
# n_species <- 6
# toy_data_set <- toy_data(n_species = n_species,
#                          times = seq(from = 0, to = 100, by = 0.1),
#                          inits = rep(0, n_species))
# 
# plot_dynamics(toy_data_set, facet = TRUE)
# 
# 
# data_collection <- get_species(data_collection, toy_data_set, 
#                                species = c(5))


# Save data
# save(data_collection, file = "lotka_volterra/toy_data_set")
load(file = "lotka_volterra/toy_data_set")

# data_collection[, 2:ncol(data_collection)] <- apply(data_collection[, 2:ncol(data_collection)],
#       2, FUN = function(i) { i/sd(i) }
#       )



## Stan Models ************************************* ####

HMM <- stan_model(file = "stan_models/HMM_gaussian.stan")
OUP_mixture <- stan_model(file = "mixture_of_experts/fit_oup_mixture.stan") 
Euler_Maruyama <- stan_model(file = "euler_maruyama/euler_maruyama.stan")

chosen_one <- "S18"

## HMM ********************************************* ####
HMM_samples <- sampling(HMM, data = list(T = length(data_collection$time),
                          K = 2, 
                          y = data_collection[, chosen_one]), 
                        iter = iter,
                        chains = chains)

hmm_parameters <- c("A", "mu", "sigma")
hmm_res <- summary(HMM_samples)$summary[, c("50%", "25%", "75%")]

hmm_res <- hmm_res[grep(paste0(hmm_parameters, collapse = "|"), rownames(hmm_res)), ] %>%
  as.data.frame() %>% 
  set_colnames(c("mode", "lower25", "upper75"))
  

hmm_map <- c(mu1 = hmm_res["mu[1]", "mode"],
             mu2 = hmm_res["mu[2]", "mode"],
             sigma1_low = hmm_res["mu[1]", "mode"] - hmm_res["sigma[1]", "mode"],
             sigma1_high = hmm_res["mu[1]", "mode"] + hmm_res["sigma[1]", "mode"],
             sigma2_low = hmm_res["mu[2]", "mode"] - hmm_res["sigma[2]", "mode"],
             sigma2_high = hmm_res["mu[2]", "mode"] + hmm_res["sigma[2]", "mode"])


hmm_plot <- ggplot(data_collection) +
  geom_ribbon(aes(x = time,
                  ymin = hmm_map["sigma1_low"],
                  ymax = hmm_map["sigma1_high"]), 
              fill = "grey") +
  geom_ribbon(aes(x = time,
                  ymin = hmm_map["sigma2_low"],
                  ymax = hmm_map["sigma2_high"]),
              fill = "grey") +
geom_line(data = data_collection, aes(x = time, y = get(chosen_one))) +
geom_hline(yintercept = c(hmm_map["mu1"], hmm_map["mu2"]),
           linetype = "dashed") +
  labs(x = "Time", y = "", subtitle = "MAP estimates for latent means +- sd", title = "HMM")



  


## OUP mixture
## Non-linear SDE
## Transformation   


## OUP mixture ************************************* ####


OUP_mixture_samples <- sampling(OUP_mixture, data = list(N = length(data_collection$time),
                                                         x = data_collection$time, 
                                                         y = data_collection[, chosen_one]), 
                                iter = iter,
                                chains = chains)


oup_mixture_parameters <- c("mu", "sigma", "lambda", "mix_weight")
oup_mixture_res <- summary(OUP_mixture_samples)$summary[, c("50%", "25%", "75%")]

oup_mixture_res <- oup_mixture_res[grep(paste0(oup_mixture_parameters, collapse = "|"), rownames(oup_mixture_res)), ] %>%
  as.data.frame() %>% 
  set_colnames(c("mode", "lower25", "upper75"))

oup_mixture_weight_res <- oup_mixture_res[grep("mix_weight", rownames(oup_mixture_res)), ] %>%
  as.data.frame() %>% 
  set_colnames(c("mode", "lower25", "upper75"))


oup_mixture_map <- c(mu1 = oup_mixture_res["mu[1]", "mode"],
                     mu2 = oup_mixture_res["mu[2]", "mode"],
                     sigma1_low = oup_mixture_res["mu[1]", "mode"] - oup_mixture_res["sigma[1]", "mode"],
                     sigma1_high = oup_mixture_res["mu[1]", "mode"] + oup_mixture_res["sigma[1]", "mode"],
                     sigma2_low = oup_mixture_res["mu[2]", "mode"] - oup_mixture_res["sigma[2]", "mode"],
                     sigma2_high = oup_mixture_res["mu[2]", "mode"] + oup_mixture_res["sigma[2]", "mode"])


oup_mixture_plot <- ggplot() +
  geom_ribbon(data = data_collection, aes(x = time,
                  ymin = oup_mixture_map["sigma1_low"],
                  ymax = oup_mixture_map["sigma1_high"]), 
              fill = "grey") +
  geom_ribbon(data = data_collection, aes(x = time,
                  ymin = oup_mixture_map["sigma2_low"],
                  ymax = oup_mixture_map["sigma2_high"]),
              fill = "grey") +
  geom_line(data = data_collection, aes(x = time, y = get(chosen_one))) +
  geom_hline(yintercept = c(oup_mixture_map["mu1"], oup_mixture_map["mu2"]),
             linetype = "dashed") +
    geom_ribbon(data = cbind(time = data_collection$time, oup_mixture_weight_res),
              aes(x = time, ymin = lower25*oup_mixture_map["mu2"], ymax = upper75*oup_mixture_map["mu2"]), fill = "firebrick") +
    geom_line(data = cbind(time = data_collection$time, oup_mixture_weight_res),
              aes(x = time, y = mode*oup_mixture_map["mu2"]), color = "red") +
  labs(x = "Time", y = "", subtitle = "MAP estimates for latent means +- marginal sd; red = mixture weight", title = "OUP Mixture")

oup_mixture_plot

## Euler-Maruyama ********************************** ####

thin_df <- data_collection[(data_collection$time %%1 == 0), ]

Euler_Maruyama_samples <- sampling(Euler_Maruyama,
                                   data = list(N = length(thin_df$time),
                                                         x = thin_df$time, 
                                                         y = thin_df[, chosen_one]), 
                                   iter = iter,
                                   chains = chains)

euler_parameters <- c("alpha", "beta", "lambda", "r", "epsilon")

euler_res <- summary(Euler_Maruyama_samples)$summary[grep(paste0(euler_parameters, collapse = "|"), rownames(summary(Euler_Maruyama_samples)$summary)), c("50%", "25%", "75%")] %>% 
  as.data.frame() %>% 
  set_colnames(c("mode", "lower25", "upper75"))
  


euler_plot <- ggplot() +
  geom_line(data = data_collection, aes(x = time, y = get(chosen_one))) +
  geom_hline(yintercept = c(euler_res[1:3, "mode"]),
             linetype = "dashed")

euler_plot

## DUMP ******************************************** ####

## fourth order polynomials
xseq <- seq(from = -1.5, to = 1.5, length.out = 100)

lapply(1:10/5, function(c) {
  
  x <- xseq
  h <- 1
  -.25*x^2*h^4 + .5*c^2*x^4
  
}) %>% do.call(cbind, .) %>%
  cbind(x = xseq, .) %>% 
  as.data.frame() %>% 
  set_colnames(c("x", paste0("V", 1:5))) %>% 
  melt(id.vars = "x") %>% 
  ggplot(aes(x = x, y = value, color =variable)) +
  geom_line()
