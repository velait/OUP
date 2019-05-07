#### Cusp on david data

## Data *************************************************** ####

load("/Users/villelaitinen/Desktop/PhD/early_warning_signals/data-David2014/David_phyloseq.Rdata")

seqA.genus <- seqA %>% aggregate_taxa("Genus")

seqA.genus.clr <- seqA.genus %>% transform("clr")

seqA.genus.abn <- seqA.genus %>% abundances() %>% t


N_genus <- 10
david_seq <- 1:100
seqA.top <- seqA.genus.abn[, top_taxa(seqA.genus, N_genus)] %>% as.data.frame() %>% mutate(time = 1:nrow(seqA.top))
  

seqA.clr.top <- (seqA.genus.clr %>% abundances() %>% t)[, top_taxa(seqA.genus, N_genus)] %>% as.data.frame() %>% mutate(time = 1:nrow(seqA.top))

seqA.top_plot <- seqA.clr.top %>%
  melt(id.vars = "time") %>% 
  ggplot(aes(x = time, y = value, color = variable)) + 
  geom_line() +
  labs(title = "Top 10 genera", subtitle = "David et al: Subject A")


seqA.top.first <- seqA.top[david_seq, ]
seqA.clr.top.first <- seqA.clr.top[david_seq, ]


# Read counts
seqA.top.first_read_counts <- seqA.top.first %>% 
  select(-time) %>% 
  rowSums()

david_stan_data <- list(N_obs = nrow(seqA.top.first),
                        N_OTUs = ncol(seqA.top.first) - 1,
                        x = seqA.top.first$time,
                        y = t(select(seqA.top.first, -time)))

## Stan *************************************************** ####

cusp_plus_obsevations <- stan_model("cusp/cusp_plus_observation.stan")

david_cusp_samples <- sampling(cusp_plus_obsevations, david_stan_data,
                                          iter = 1000, 
                                          chains = 1, 
                                          control = list(adapt_delta = 0.8))


## Results ************************************************ ####


# get results
results_df <- lapply(1:4, function(x) get_stan_results(david_cusp_samples, paste0("theta\\[.+,", x,"\\]"))) %>% do.call(rbind, .) %>% 
  cbind(parameter = rep(c("alpha", "beta", "lambda", "epsilon"), each = N),
        index = factor(rep(1:N, 4)))


# Map estimates
david_map <- lapply(1:10, function(i) {
  
  rr <- 1
  aa <- results_df %>%
    filter(parameter == "alpha", index == i) %>% 
    pull(mode)
  
  bb <- results_df %>%
    filter(parameter == "beta", index == i) %>% 
    pull(mode)
  
  ee <- results_df %>%
    filter(parameter == "epsilon", index == i) %>% 
    pull(mode)
  
  ll <- results_df %>%
    filter(parameter == "lambda", index == i) %>% 
    pull(mode)
  
  d <- cc_density(x, rr, aa, bb, ll, ee)
  
  data.frame(x, d, i) %>% 
    set_colnames(c("x", "y", "index"))
  
}) %>% do.call(rbind, .)


# Compositional data densities

first_cs <- david_stan_data$y %>% colSums()

compositional_first <- lapply(1:ncol(david_stan_data$y), function(i) {
  
  david_stan_data$y[, i]/first_cs[i]
  
}) %>% do.call(cbind, .) %>% as.data.frame()


compositional_densities <- compositional_first %>% t %>% melt %>% mutate(index = rep(1:10, each = 100))

# Plot
p <- ggplot() +
  geom_line(data = david_map, aes(x = x, y = y), color = "red") +
  facet_wrap(~index) +
  labs(title = "David MAP estimates")
  

q <- ggplot() +
  stat_density(data = compositional_densities, aes(x = value), geom = "line") +
  facet_wrap(~index) +
  labs(title = "Corresponding compositional densities")


bug_enumeration <- c(top_taxa(seqA.genus, N_genus)) %>% `names<-`(1:N_genus %>% as.character())


## Get samples from stan object, simulate and plot ******** ####


n_sim <- 20

cusp_parameters <- c("alpha", "beta", "lambda", "epsilon")

simulated_densities <- lapply(1, function(i) {
  
  # number of samples
  n_samples <- david_cusp_samples@stan_args[[1]]$iter - david_cusp_samples@stan_args[[1]]$warmup
  
  # samples indices
  sample_ind <- sample(1:n_samples, n_sim)
  
  par_samples <- rstan::extract(david_cusp_samples, "theta")[[1]][sample_ind, , ]
  
  
  
  
  # simulate series
  series_sim <- lapply(1:n_sim, function(s) {
    lapply(1:N_genus, function(j) {
      shoji_generator(y0 = rnorm(1, 0, 1), times = seq(from=david_seq[1], to=david_seq[length(david_seq)], by = 0.01),
                      r = 1,
                      alpha = par_samples[1, j, 1],
                      beta = par_samples[1, j, 2],
                      lambda = par_samples[1, j, 3],
                      epsilon = par_samples[1, j, 4],
                      seed = sample(1:1000, 1))
    }) %>%
      do.call(cbind, .) %>%
      as.data.frame() %>%
      cbind(time = seq(from=david_seq[1], to=david_seq[length(david_seq)], by = 0.01)) %>%
      set_colnames(c(paste0("y", 1:N_genus), "x")) %>% 
      mutate(sample = s)
  }) %>% 
    do.call(rbind, .)
    

  
  # thin: dt = 1
  series_sim <- thin(series_sim)
  series_sim_no_time <- series_sim %>% select(-c(x, sample))

  
  # Softmax at each time step
  series_sim_softmax <- lapply(1:nrow(series_sim), function(r) {
    softMax(series_sim_no_time[r, ])
  }) %>%
    do.call(rbind, .) %>% 
    cbind(., series_sim %>%
            select(c(x, sample)))
    
  
  series_sim_softmax_no_time <- series_sim_softmax %>% select(-c(x, sample))
  
  
  # Multinomial samping; read counts from data
  david_reads <- rep(seqA.top.first_read_counts, each = length(david_seq))
  series_sim_multin <- lapply(1:nrow(series_sim_softmax_no_time), function(r) {
    rmultinom(n = 1, size = david_reads[r], prob = series_sim_softmax_no_time[r, ]) %>% 
      t
  }) %>% 
    do.call(rbind, .) %>% 
    cbind(., series_sim %>%
            select(c(x, sample)))
  
  
  rm(david_reads, series_sim_softmax_no_time, series_sim_no_time)
  
  return(list(cusp = series_sim, softmax = series_sim_softmax, multin = series_sim_multin))
  
})


# Transform to compositional
simulated_multin <- simulated_densities[[1]][["multin"]] %>% 
  select(-c(x, sample))

simulated_multin_compositional <- lapply(1:nrow(simulated_multin), function(i) {
  
  simulated_multin[i, ]/rowSums(simulated_multin)[i]
  
}) %>%
  do.call(rbind, .) %>% 
  cbind(simulated_densities[[1]][["multin"]][, c("x", "sample")])

# Plot
ggplot() + 
  stat_density(data =  simulated_multin_compositional %>% 
                 melt(id.vars = c("x", "sample")),
               aes(x = value, group = sample),
               geom = "line", position = "identity", alpha = 0.5, color = "red") +
  stat_density(data = compositional_densities %>% mutate(sample = index), 
            aes(x = value), geom = "line") +
  facet_wrap(~variable, scales = "free") +
  labs(subtitle = "David et al data and posterior samples")



## CLR abundance; no observation model ******************** ####



david_clr_stan_data <- list(N_obs = nrow(seqA.clr.top.first),
                        N_OTUs = ncol(seqA.clr.top.first) - 1,
                        x = seqA.top.first$time,
                        y = t(select(seqA.clr.top.first, -time)))

## Stan *************************************************** ####

cusp_clr <- stan_model("cusp/shoji_cusp.stan")

david_cusp_clr_samples <- sampling(cusp_clr,
                               david_clr_stan_data,
                               iter = 1000, 
                               chains = 1, 
                               control = list(adapt_delta = 0.8))


## Results ************************************************ ####


# get results
clr_results_df <- lapply(1:4, function(x) get_stan_results(david_cusp_clr_samples
, paste0("theta\\[.+,", x,"\\]"))) %>% do.call(rbind, .) %>% 
  cbind(parameter = rep(c("alpha", "beta", "lambda", "epsilon"), each = N),
        index = factor(rep(1:N, 4)))


# Map estimates
david_clr_map <- lapply(1:10, function(i) {
  
  rr <- 1
  aa <- clr_results_df %>%
    filter(parameter == "alpha", index == i) %>% 
    pull(mode)
  
  bb <- clr_results_df %>%
    filter(parameter == "beta", index == i) %>% 
    pull(mode)
  
  ee <- clr_results_df %>%
    filter(parameter == "epsilon", index == i) %>% 
    pull(mode)
  
  ll <- clr_results_df %>%
    filter(parameter == "lambda", index == i) %>% 
    pull(mode)
  
  x <- seq(from = 0, to = 15,  length.out = 100)
  
  d <- cc_density(x, rr, aa, bb, ll, ee)
  
  data.frame(x, d, i) %>% 
    set_colnames(c("x", "y", "index"))
  
}) %>% do.call(rbind, .)


# Plot
p <- ggplot() +
  geom_line(data = david_clr_map, aes(x = x, y = y), color = "red") +
  facet_wrap(~index) +
  labs(title = "David (CLR)  MAP estimates")


q <- ggplot() +
  stat_density(data = seqA.clr.top.first %>% melt(id.vars = "time"), aes(x = value), geom = "line") +
  facet_wrap(~variable) +
  labs(title = "Corresponding compositional densities")


bug_enumeration <- c(top_taxa(seqA.genus, N_genus)) %>% `names<-`(1:N_genus %>% as.character())

ggplot() + 
  stat_density(data = seqA.clr.top.first %>%
                 set_colnames(c(1:10 %>% as.character, "time")) %>% 
                 melt(id.vars = "time") %>%
                 mutate(index = variable),
               aes(x = value), geom = "line") +
  geom_line(data = david_clr_map, aes(x = x, y = y), color = "red") +
  facet_wrap(~index, labeller = labeller(index = bug_enumeration)) +
  labs(title = "David et al data (CLR) and posterior samples",
       subtitle = "Red = MAP, Black = data; 100 time points")


## Get samples from stan object, simulate and plot ******** ####


n_sim <- 10

cusp_parameters <- c("alpha", "beta", "lambda", "epsilon")

simulated_densities_clr <- lapply(1, function(i) {
  
  # number of samples
  n_samples <- david_cusp_clr_samples@stan_args[[1]]$iter - david_cusp_clr_samples@stan_args[[1]]$warmup
  
  # samples indices
  sample_ind <- sample(1:n_samples, n_sim)
  
  par_samples <- rstan::extract(david_cusp_clr_samples, "theta")[[1]][sample_ind, , ]
  
  
  
  # simulate series
  ticks <- 0.1
  series_sim <- lapply(1:n_sim, function(s) {
    lapply(1:N_genus, function(j) {
      shoji_generator(y0 = rnorm(1, 5, 1), times = seq(from=david_seq[1], to=david_seq[length(david_seq)], by = ticks),
                      r = 1,
                      alpha = par_samples[1, j, 1],
                      beta = par_samples[1, j, 2],
                      lambda = par_samples[1, j, 3],
                      epsilon = par_samples[1, j, 4],
                      seed = sample(1:1000, 1))
    }) %>%
      do.call(cbind, .) %>%
      as.data.frame() %>%
      cbind(time = seq(from=david_seq[1], to=david_seq[length(david_seq)], by = ticks)) %>%
      set_colnames(c(paste0("y", 1:N_genus), "x")) %>% 
      mutate(sample = s)
  }) %>% 
    do.call(rbind, .)
  
  
  
  # thin: dt = 1
  series_sim <- thin(series_sim)
  
  return(series_sim)
  
})

# Plot
ggplot() +
  stat_density(data = simulated_densities_clr[[1]] %>%
                 set_colnames(c(1:10 %>% as.character, "x", "sample")) %>%
                 melt(id.vars = c("x", "sample")),
               aes(x = value, color = sample),
               geom = "line", position = "identity", color = "grey") +
  stat_density(data =  seqA.clr.top.first %>%
                 set_colnames(c(1:10 %>% as.character, "time")) %>%
                 melt(id.vars = "time"),
               aes(x = value, group = variable),
               geom = "line", position = "identity") +
  facet_wrap(~variable, scales = "free", labeller = labeller(variable = bug_enumeration)) +
  labs(subtitle = "David et al data (CLR). Grey = posterior samples --> averaged stationary density") 

