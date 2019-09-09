# Hierarchical Cusp model


## Parameters  **************** ####

seed <- 11235
set.seed(seed)


N_obs <- 10
N_OTUs <- 50

times <- seq(from = 1, to = N_obs, by = 1)


cusp_hp <- list(alpha = c(0, 1),
                beta = c(0, 1),
                lambda = c(0, 1),
                epsilon = c(0, 1),
                # r = c(1, 1)
                r = c(4, 4)
                )

cusp_parameters <- list(alpha = rnorm(N_OTUs, cusp_hp[["alpha"]][1], cusp_hp[["alpha"]][2]),
                        beta = rnorm(N_OTUs, cusp_hp[["beta"]][1], cusp_hp[["beta"]][2]),
                        lambda = rnorm(N_OTUs, cusp_hp[["lambda"]][1], cusp_hp[["lambda"]][2]),
                        epsilon = rnorm_trunc(N_OTUs, cusp_hp[["epsilon"]][1], cusp_hp[["epsilon"]][2]),
                        # r = rnorm_trunc(N_OTUs, cusp_hp[["r"]][1], cusp_hp[["r"]][2])
                        # r = rgamma(N_OTUs, cusp_hp[["r"]][1], cusp_hp[["r"]][2])
                        r = rep(1, N_OTUs)
                        )


## Data *********************** ####

hierarchical_cusp_obs <- lapply(1:N_OTUs, function(otu) {
  
  alpha <- cusp_parameters[["alpha"]][otu]
  beta <- cusp_parameters[["beta"]][otu]
  lambda <- cusp_parameters[["lambda"]][otu]
  epsilon <- cusp_parameters[["epsilon"]][otu]
  r <- cusp_parameters[["r"]][otu]
  
  
  x <- shoji_generator(y0 = rnorm(1, lambda, epsilon),
                  times = seq(from = times[1], to = times[length(times)], by = 0.1),
                  r = r,
                  alpha = alpha,
                  beta = beta,
                  lambda = lambda,
                  epsilon = epsilon,
                  seed = NULL)

  
  return(x)
  
}) %>%
  do.call(cbind, .) %>% 
  cbind(x = seq(from = times[1], to = times[length(times)], by = 0.1)) %>% 
  thin() %>% 
  as.data.frame()


hierarchical_cusp_stan_data <- list(N_OTUs = N_OTUs, 
                                    N_obs, 
                                    x = hierarchical_cusp_obs$x,
                                    y = hierarchical_cusp_obs %>% select(-x) %>% t)

## Stan  ********************** ####

## Models
# Regular
# Hierarchical
hierarchical_cusp_model <- stan_model("cusp/hierarchical_cusp/hierarchical_cusp_shoji.stan")
cusp_model <- stan_model("cusp/shoji_cusp_with_r.stan")
tied_model <- stan_model("cusp/hierarchical_cusp/cusp_tied_shoji.stan")

## Samples
# Hierarchical
hierarchical_cusp_samples <- sampling(hierarchical_cusp_model, 
                                      hierarchical_cusp_stan_data, 
                                      iter = 1000, 
                                      chains = 1)

# Regular
cusp_samples <- sampling(cusp_model, 
                         hierarchical_cusp_stan_data, 
                         iter = 1000, 
                         chains = 1)


# Tied
tied_samples <- sampling(tied_model, 
                         hierarchical_cusp_stan_data, 
                         iter = 1000, 
                         chains = 1)

## Results ******************** ####

N_samples <- 30


## Parameter estimates ********

# Hierarchical and regular summary statistics
par_res <- lapply(names(cusp_parameters), function(par) {
  
  # Hierarchical
  hier <- get_stan_results(hierarchical_cusp_samples, paste0(par, "\\["))
  
  hier_df <- hier %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "parameter") %>% 
    mutate(true = cusp_parameters[[par]], 
           index =  str_extract_all(parameter, "[0-9]+") %>% as.numeric(), 
           parameter = par,
           model = "hierarchical")
    
  
  
  # Regular
  reg <- get_stan_results(cusp_samples, paste0(par, "\\["))
  reg_df <- reg %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "parameter") %>% 
    mutate(true = cusp_parameters[[par]], 
           index =  str_extract_all(parameter, "[0-9]+")%>% as.numeric(), 
           parameter = par,
           model = "regular")
  
  
  
  df <- rbind(hier_df, reg_df) %>% 
    arrange(true)
  
  
  return(df)
  
}) %>% 
  do.call(rbind, .)

# Tied parameters summary statistics
tied_res <- lapply(names(cusp_parameters), function(par) {
  
  res <- get_stan_results(tied_samples, paste0(par))
  
  res_df <- res %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "parameter") %>% 
    mutate(parameter = par,
           model = "tied")
  
}) %>% do.call(rbind, .)

# Samples from individual parameter posteriors, hierarchical and regular
sample_res <- lapply(1, function(i) {
  
  sample_ind <- sample(1:(rstan::extract(hierarchical_cusp_samples, par)[[1]] %>% dim)[1], N_samples)
  
  hier_df <- lapply(names(cusp_parameters), function(par) {
    
    rstan::extract(hierarchical_cusp_samples, par)[[1]][sample_ind, ] %>%
      as.vector() %>% 
      cbind(rep(1:N_samples, N_OTUs), 
            rep(1:N_OTUs, each = N_samples)) %>%
      as.data.frame() %>% 
      mutate(parameter = par) %>% 
      set_colnames(c("value", "sample", "index", "parameter"))
  }) %>% do.call(rbind, .) %>% mutate(model = "hierarchical")
  
  
  
  reg_df <- lapply(names(cusp_parameters), function(par) {
    
    rstan::extract(cusp_samples, par)[[1]][sample_ind, ] %>%
      as.vector() %>% 
      cbind(rep(1:N_samples, N_OTUs), 
            rep(1:N_OTUs, each = N_samples)) %>%
      as.data.frame() %>% 
      mutate(parameter = par) %>% 
      set_colnames(c("value", "sample", "index", "parameter"))
  }) %>% do.call(rbind, .)  %>% mutate(model = "regular")
  
  
  
  return(rbind(hier_df, reg_df))
  
}) %>% do.call(rbind, .)
  
  
# Plot

par_plots <- lapply(names(cusp_parameters), function(par) {
  
  mydf <- par_res %>%
    filter(parameter == par) %>% 
    arrange(true) %>% 
    mutate(dummy = rep(1:N_OTUs, each = 2))
  
  p <- ggplot() +
    geom_line(data = mydf, aes(x = dummy, y = true)) +
    geom_hline(data = tied_res %>% filter(parameter == par), aes(yintercept = mode), linetype = "dashed") + 
    geom_errorbar(data = mydf, aes(x = dummy, ymin = lower2.5, ymax = upper97.5, color = model), position = "dodge") +
    labs(title = par) +
    scale_color_tron()
  
  p
})

plot_grid(par_plots[[1]], 
          par_plots[[2]], 
          par_plots[[3]], 
          par_plots[[4]], 
          par_plots[[5]])

## Population distributions *******

hyperpar_res <- lapply(names(cusp_parameters), function(par) {
  
  res <- get_stan_results(hierarchical_cusp_samples, paste0(par, "_")) %>% 
    rownames_to_column(var = "parameter") %>% 
    mutate(index =  str_extract_all(parameter, "[0-9]+") %>% as.numeric(), 
           parameter = par)
  
  # Add true values
  res <- res %>% 
    mutate(true = cusp_hp[[par]], 
           type = "MAP", 
           sample = NA) %>% 
    select(parameter, mode, index, sample, true, type) %>% 
    set_colnames(c("parameter", "value", "index", "sample",  "true", "type"))
  
  # Samples of hyperparameters for the paramter par
  hp_samples <- cbind(rstan::extract(hierarchical_cusp_samples, paste0(par, "_1"))[[1]],
  rstan::extract(hierarchical_cusp_samples, paste0(par, "_2"))[[1]])
  
  
  ss <- hp_samples[sample(1:nrow(hp_samples), N_samples), ] %>%
    as.vector() %>% 
    cbind(rep(1:2, each = N_samples), 
          rep(1:N_samples, 2)) %>% 
    as.data.frame() %>% 
    mutate(parameter = par) %>% 
    set_colnames(c("value", "index", "sample", "parameter")) %>%  
    mutate(true = NA, type = "sample") %>% 
    select(parameter, value, index, sample, true, type)
  
  
  
  
  return(rbind(res, ss))
  
  
}) %>% do.call(rbind, .)


population_sample_densities <- lapply(1:N_samples, function(i) {
  
  # Get sample
  df <- hyperpar_res %>% 
    filter(sample == i)
  
  # x_grid
  x <- seq(from = -5, to = 5, by = 0.1)
  
  # distribution values for each hyper paramter
  value_df <- lapply(names(cusp_parameters), function(par) {
    
    vals <- df %>% filter(parameter == par) %>% pull(value)
    
    data.frame(y = dnorm(x, vals[1], vals[2]), x = x, parameter = par)
    
  }) %>% 
    do.call(rbind, .) %>% 
    cbind(sample = i)
  
  
  
  return(value_df)
}) %>% 
  do.call(rbind, .)

population_map_densities <- lapply(1, function(i) {
  
  df <- hyperpar_res %>% 
    filter(type == "MAP")
  
  # x grid
  x <- seq(from = -5, to = 5, by = 0.1)
  
  
  
  # distribution values for each hyper paramter
  value_df <- lapply(names(cusp_parameters), function(par) {
    
    vals <- df %>% filter(parameter == par) %>% pull(value)
    
    data.frame(y = dnorm(x, vals[1], vals[2]), x = x, parameter = par)
    
  }) %>% 
    do.call(rbind, .)
  
  
  return(value_df)
})[[1]]


population_true_densities <- lapply(1, function(i) {
  
  df <- hyperpar_res %>% 
    filter(type == "MAP")
  
  # x grid
  x <- seq(from = -5, to = 5, by = 0.1)
  
  
  
  # distribution values for each hyper paramter
  value_df <- lapply(names(cusp_parameters), function(par) {
    
    vals <- df %>% filter(parameter == par) %>% pull(true)
    
    data.frame(y = dnorm(x, vals[1], vals[2]), x = x, parameter = par)
    
  }) %>% 
    do.call(rbind, .)
  
  
  return(value_df)
})[[1]]


# Plot
ggplot() +
  geom_line(data = population_sample_densities, aes(x = x, y = y, group = sample), color = "salmon", alpha = 0.1) +
  geom_line(data = population_map_densities, aes(x = x, y = y), color = "red") + 
  geom_line(data = population_true_densities, aes(x = x, y = y)) + 
  facet_wrap(~parameter, scales = "free") + 
  labs(title = "Population distributions", subtitle = "Black = true")

## Simulated densities *************

#hierarchical
hier_simulated_series <- lapply(1:N_OTUs, function(i) {
  
  df <- sample_res %>% 
    filter(index == i, model == "hierarchical")
  
  
  sim_df <- lapply(1:N_samples, function(s) {
    
    s_df <- df %>% 
      filter(sample == s)
    
    as <- s_df %>% filter(parameter == "alpha") %>% pull(value)
    bs <- s_df %>% filter(parameter == "beta") %>% pull(value)
    ls <- s_df %>% filter(parameter == "lambda") %>% pull(value)
    es <- s_df %>% filter(parameter == "epsilon") %>% pull(value)
    rs <- s_df %>% filter(parameter == "r") %>% pull(value)
    
    x <- seq(from = 0, to = 100, by = 0.1)
    
    series_df <- shoji_generator(y0 = rnorm(1, ls, es), times = x,
                    r = rs, alpha = as, beta = bs, lambda = ls, epsilon = es, seed = NULL) %>% 
      data.frame(y = ., x = x) %>% 
      thin() %>% 
      cbind(sample = s)
    
  }) %>% do.call(rbind, .)
  
  
  # Add otu indentifier
  sim_df <- sim_df %>% 
    mutate(index = i)
  
  return(sim_df)
}) %>% 
  do.call(rbind, .)

hier_map_simulated_series <- lapply(1:N_OTUs, function(i) {
  
  df <- par_res %>% 
    filter(index == i, model == "hierarchical")
  
  
  sim_df <- lapply(1, function(s) {
    
    as <- df %>% filter(parameter == "alpha") %>% pull(mode)
    bs <- df %>% filter(parameter == "beta") %>% pull(mode)
    ls <- df %>% filter(parameter == "lambda") %>% pull(mode)
    es <- df %>% filter(parameter == "epsilon") %>% pull(mode)
    rs <- df %>% filter(parameter == "r") %>% pull(mode)
    
    x <- seq(from = 0, to = 100, by = 0.1)
    
    series_df <- shoji_generator(y0 = rnorm(1, ls, es), times = x,
                                 r = rs, alpha = as, beta = bs, lambda = ls, epsilon = es, seed = NULL) %>% 
      data.frame(y = ., x = x) %>% 
      thin()
    
  }) %>% do.call(rbind, .)
  
  
  # Add otu indentifier
  sim_df <- sim_df %>% 
    mutate(index = i)
  
  return(sim_df)
}) %>% 
  do.call(rbind, .)

#regular
reg_simulated_series <- lapply(1:N_OTUs, function(i) {
  
  df <- sample_res %>% 
    filter(index == i, model == "regular")
  
  
  sim_df <- lapply(1:N_samples, function(s) {
    
    s_df <- df %>% 
      filter(sample == s)
    
    as <- s_df %>% filter(parameter == "alpha") %>% pull(value)
    bs <- s_df %>% filter(parameter == "beta") %>% pull(value)
    ls <- s_df %>% filter(parameter == "lambda") %>% pull(value)
    es <- s_df %>% filter(parameter == "epsilon") %>% pull(value)
    rs <- s_df %>% filter(parameter == "r") %>% pull(value)
    
    x <- seq(from = 0, to = 100, by = 0.1)
    
    series_df <- shoji_generator(y0 = rnorm(1, ls, es), times = x,
                                 r = rs, alpha = as, beta = bs, lambda = ls, epsilon = es, seed = NULL) %>% 
      data.frame(y = ., x = x) %>% 
      thin() %>% 
      cbind(sample = s)
    
  }) %>% do.call(rbind, .)
  
  
  # Add otu indentifier
  sim_df <- sim_df %>% 
    mutate(index = i)
  
  return(sim_df)
}) %>% 
  do.call(rbind, .)

reg_map_simulated_series <- lapply(1:N_OTUs, function(i) {
  
  df <- par_res %>% 
    filter(index == i, model == "regular")
  
  
  sim_df <- lapply(1, function(s) {
    
    as <- df %>% filter(parameter == "alpha") %>% pull(mode)
    bs <- df %>% filter(parameter == "beta") %>% pull(mode)
    ls <- df %>% filter(parameter == "lambda") %>% pull(mode)
    es <- df %>% filter(parameter == "epsilon") %>% pull(mode)
    rs <- df %>% filter(parameter == "r") %>% pull(mode)
    
    x <- seq(from = 0, to = 100, by = 0.1)
    
    series_df <- shoji_generator(y0 = rnorm(1, ls, es), times = x,
                                 r = rs, alpha = as, beta = bs, lambda = ls, epsilon = es, seed = NULL) %>% 
      data.frame(y = ., x = x) %>% 
      thin()
    
  }) %>% do.call(rbind, .)
  
  
  # Add otu indentifier
  sim_df <- sim_df %>% 
    mutate(index = i)
  
  return(sim_df)
}) %>% 
  do.call(rbind, .)

true_densities <- lapply(1:N_OTUs, function(i) {
  
  
  
  
  as <- cusp_parameters[["alpha"]][i]
  bs <- cusp_parameters[["beta"]][i]
  ls <- cusp_parameters[["lambda"]][i]
  es <- cusp_parameters[["epsilon"]][i]
  rs <- cusp_parameters[["r"]][i]
  
  x <- seq(from = -5, to = 5, by = 0.1)
  
  series_df <- cc_density(x = x, rr = rs, aa = as, bb =  bs, ll = ls, ee = es) %>% 
    cbind(index = i) %>%
    as.data.frame() %>% 
    mutate(x = x) %>% 
    set_colnames(c("y", "index", "x"))

  return(series_df)
}) %>% 
  do.call(rbind, .)



ggplot() + 
  geom_line(data = true_densities, aes(x = x, y = y), size = 1) +
  stat_density(data = reg_simulated_series, aes(x = y, group = sample), alpha = 0.1, color = "royalblue", geom = "line", position = "identity") +
  stat_density(data = reg_map_simulated_series, aes(x = y), color = "blue", geom = "line", position = "identity") +
  stat_density(data = hier_simulated_series, aes(x = y, group = sample), alpha = 0.1, color = "salmon", geom = "line", position = "identity") +
  stat_density(data = hier_map_simulated_series, aes(x = y), color = "red", geom = "line", position = "identity") +
  stat_density(data = hierarchical_cusp_obs %>%
                 melt(id.vars = "x") %>% 
                 mutate(index = variable %>%
                          gsub("V", "", .),
                        y = value),
               aes(x = y),
               color = "black", geom = "line", position = "identity", linetype = "dashed") + 
  facet_wrap(~index, scales = "free", ncol = 5) +
  labs(subtitle = "Red = hiearchical; Blue = regular; Dashed = data; Solid = true")






