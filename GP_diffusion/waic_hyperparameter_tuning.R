## WAIC hyperparameter tuning for a taxon

hitchip_otu <- "Akkermansia"


length_scale_par <- 2
variance_par <- .25
epsilon <- .1


x_pred <- range_seq(range(hitchip_data$abundances[, hitchip_otu])[1:2], length_out = 1)
mean_c <- mean(hitchip_data$abundances[, hitchip_otu])

## One dimensional
data <- hitchip_data$delta_df %>% 
  filter(otu == hitchip_otu) %>% 
  filter(dx != 0) 

## Grid search values
rho_grid <- 1:10/2
sigma_grid <- 1:10/10
epsilon_grid <- 1:10/10

hyperparameter_grid <- expand.grid(rho_grid, sigma_grid, epsilon_grid) %>% 
  set_colnames(c("rho", "sigma", "epsilon"))


waic_taxa <- test_taxa
waic_results <- data.frame(NA, NA, NA, NA, NA, NA) %>% 
  set_colnames(c("rho", "sigma", "epsilon" ,"WAIC", "target", "OTU"))

for(taxon in waic_taxa) {
  
  ## Drift
  drift_tuning <- lapply(1:nrow(hyperparameter_grid), function(i) {
    print(i)
    hp <- hyperparameter_grid[i, ]
    length_scale_par <- hp$rho
    variance_par <- hp$sigma
    obs_noise <- hp$epsilon
    
    drift_imputation <- GP_imputation(X = data %>% 
                                        pull(x), 
                                      Y = data %>%
                                        pull(unit_dx),
                                      X_pred = x_pred,
                                      rho = length_scale_par,
                                      sigma = variance_par,
                                      epsilon = obs_noise, 
                                      mean_c = mean_c,
                                      add_linear = FALSE, 
                                      n_samples = n_samples)
    
    data.frame(hp, WAIC = drift_imputation$WAIC, target = "drift")
  }) %>%
    do.call(rbind, .)
  
  
  ## Diffusion
  diffusion_tuning <- lapply(1:nrow(hyperparameter_grid), function(i) {
    print(i)
    hp <- hyperparameter_grid[i, ]
    length_scale_par <- hp$rho
    variance_par <- hp$sigma
    obs_noise <- hp$epsilon
    
    drift_imputation <- GP_imputation(X = data %>% 
                                        pull(x), 
                                      Y = data %>%
                                        pull(dx) %>%
                                        abs,
                                      X_pred = x_pred,
                                      rho = length_scale_par,
                                      sigma = variance_par,
                                      epsilon = obs_noise, 
                                      mean_c = mean_c,
                                      add_linear = FALSE, 
                                      n_samples = n_samples)
    
    data.frame(hp, WAIC = drift_imputation$WAIC, target = "diffusion")
  }) %>%
    do.call(rbind, .)
  
  
  OTU_df <- rbind(drift_tuning, 
        diffusion_tuning) %>% 
    mutate(OTU = taxon)
  
  waic_results <- waic_results %>% 
    rbind(OTU_df)
  
  saveRDS(waic_results, file = "HitChip_tests/HitChip_log10_WAIC_results.Rdata")
  
}


# Filter optimal results per taxa
waic_results <- waic_results %>% drop_na

waic_results <- waic_results %>% 
  group_by(OTU, target) %>% 
  filter(WAIC == max(WAIC)) %>% 
  ungroup %>% 
  distinct


## WAIC tiles

waic_results %>% 
  ggplot(aes(x = ))

## Plot with optimized parameters
waic_OTUs <- waic_results$OTU %>% unique


stationary_density_plots <- lapply(waic_OTUs, function(taxon) {
  print(taxon)
  
  taxon_pars <- waic_results %>% 
    filter(OTU == taxon) %>% 
    select(rho, sigma, epsilon)
  
  drift_pars <- taxon_pars[1, ]
  diffusion_pars <- taxon_pars[2, ]
  
  p <- hitchip_stationary_density_plotter(hitchip_data, taxon, 
                                          drift_pars, diffusion_pars)
  
  return(p)
}) %>%
  set_names(waic_OTUs)



plot_grid(stationary_density_plots[[1]] + coord_cartesian(ylim = 0:2), 
          stationary_density_plots[[2]] + coord_cartesian(ylim = 0:2), 
          stationary_density_plots[[3]] + coord_cartesian(ylim = 0:2), 
          stationary_density_plots[[4]] + coord_cartesian(ylim = 0:2), 
          stationary_density_plots[[5]] + coord_cartesian(ylim = 0:2), 
          stationary_density_plots[[6]] + coord_cartesian(ylim = 0:2), 
          stationary_density_plots[[7]] + coord_cartesian(ylim = 0:2), 
          stationary_density_plots[[8]] + coord_cartesian(ylim = 0:2), 
          stationary_density_plots[[9]] + coord_cartesian(ylim = 0:2), 
          stationary_density_plots[[10]] + coord_cartesian(ylim = 0:2), 
          stationary_density_plots[[11]] + coord_cartesian(ylim = 0:2), 
          stationary_density_plots[[12]] + coord_cartesian(ylim = 0:2), 
          stationary_density_plots[[13]] + coord_cartesian(ylim = 0:2),
          stationary_density_plots[[14]] + coord_cartesian(ylim = 0:2),
          stationary_density_plots[[15]] + coord_cartesian(ylim = 0:2), 
          stationary_density_plots[[16]] + coord_cartesian(ylim = 0:2), 
          stationary_density_plots[[17]] + coord_cartesian(ylim = 0:2),
          ncol= 4)
