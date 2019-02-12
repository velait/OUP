rho_lh_grid <- 1:125/10
alpha_lh_grid <- 1:125/10
measurement_error <- sigma

# Covariance matrices
covM <- lapply(rho_lh_grid, function(r) {
  a_cov <- lapply(alpha_lh_grid, function(a) {
    
    print(paste0("rho: ", r, "/", max(rho_lh_grid), "; alpha: ", a, "/", max(alpha_lh_grid)))
    # covariance matrix
    
    covM <- matrix(NA, 101, 101)
    for(i in 1:101) {
      for(j in 1:101) {
        
        covM[i, j] <- a^2*exp(-abs(i-j)/(r^2)) + ifelse(i == j, measurement_error^2, 0)
        
      }
    }
    
    return(covM)
    
  }) %>% set_names(alpha_lh_grid %>% as.character())
  return(a_cov)
}) %>% set_names(rho_lh_grid %>% as.character())

# Prior 
print("prior")
prior <- data.frame(rho = rep(rho_lh_grid, length(alpha_lh_grid)),
                    alpha = rep(alpha_lh_grid, each = length(rho_lh_grid)),
                    value = NA)

for(alpha in alpha_lh_grid) {
  
  for(rho in rho_lh_grid) {
    
    condition <- rho == prior$rho & alpha == prior$alpha
    
    prior[condition, "value"] <- dinvgamma(rho, 4, 10)*dnorm(alpha, 0, 1)
    
  }
  
}


# Likelihoods and posteriors
manual_results <- lapply(rho_grid, function(r) {
  
  alpha_results <- lapply(alpha_grid, function(a) {
    print(paste0("rho: ", r, "/", length(rho_grid), "; alpha: ", a, "/", length(alpha_grid)))
    

    ## Likelihood ****************************** ####
    print("likelihood")
    # Fix sigma at 1
    
    rho_true <- r
    alpha_true <- a
    measurement_error <- sigma
    
    # data
    dat <- single_simulated_data[[rho_true]][[alpha_true]] %>% 
      filter((x%%1==0))
    dat <- list(N = nrow(dat), y = dat$y, x = dat$x)
    
    
    likelihood_val <- data.frame(rho = rep(rho_lh_grid, length(alpha_lh_grid)),
                                 alpha = rep(alpha_lh_grid, each = length(rho_lh_grid)),
                                 value = NA)
    
    for(alpha in alpha_lh_grid) {
      
      for(rho in rho_lh_grid) {
        
        condition <- rho == likelihood_val$rho & alpha == likelihood_val$alpha
        
        likelihood_val[condition, "value"] <- dmvnorm(x = dat$y, sigma = covM[[as.character(rho)]][[as.character(alpha)]])
        
      }
      
    }
    
    ## posterior ****************************** ####
    print("posterior")
    posterior_val <- cbind(likelihood_val, prior = prior$value) %>% 
      set_colnames(c("rho", "alpha", "likelihood", "prior")) %>% 
      mutate(posterior = likelihood*prior) %>% 
      mutate(posterior = posterior/sum(posterior)) %>%
      mutate(log_likelihood = log(likelihood)) %>% 
      melt(id.vars = c("alpha", "rho"))
    

    return(posterior_val)
    
  }) %>% set_names(alpha_grid %>% as.character)
  
  return(alpha_results)
  
}) %>% set_names(rho_grid %>% as.character)

# save(manual_results,
#      covM,
#      file = "results/manual_computation_results.Rdata")


# load(file = "results/manual_computation_results.Rdata")




## Plot ********** ####

rho_true <- 1
alpha_true <- 1
measurement_error <- 1

manual_results_distribution_plot <- lapply(rho_grid, function(r) {
  
  alpha_results <- lapply(alpha_grid, function(a) {
    
    p <- manual_results[[r]][[a]] %>% 
      ggplot() +
      geom_contour(aes(x = rho, y = alpha,  z = value, color = variable)) +
      scale_color_startrek() +
      theme_bw() +
      coord_fixed() +
      geom_abline(intercept = 0, slope = a/r, linetype = "dashed") +
      geom_point(aes(x = r, y = a)) +
      labs(x = "Rho", y = "Alpha", title = paste0("True values: rho = ", r, ", alpha = ", a))
    
    p
    
    
  }) %>% set_names(alpha_grid %>% as.character())
  
}) %>% set_names(rho_grid %>% as.character())

# Save plots as png
for(r in rho_grid) {
  for(a in alpha_grid) {
    print(a)
    
    
    png(paste0("figures/manual_computations/manual_results_distribution_plot_", r, "_", a, ".png"),
        , width = 800, height = 600)
    print(manual_results_distribution_plot[[r]][[a]])
    dev.off()
    
  }
}



## Compare true, stan and manual

# Get manual MAP
manual_MAP <- lapply(rho_grid, function(r) {
  
  lapply(alpha_grid, function(a) {
    
    df <- manual_results[[r]][[a]] %>% 
      filter(variable == "posterior")
    
    manual_MAP <-  (df[which.max(df$value), c("alpha", "rho")] %>% t)[, 1]
    
    data.frame(manual_MAP = manual_MAP,
               parameter = c("alpha", "rho"),
               true_value = c(a, r))
      
  
  })
  
})



estimate_comparison_df <- lapply(rho_grid, function(r) {
  
  lapply(alpha_grid, function(a) {
    
    df <- cbind(single_process_results[[r]][[a]][c("alpha", "rho"), c("mode", "true_value")], 
    manual_MAP[[r]][[a]][, 1:2])
    
  }) %>% do.call(rbind, .)
})%>% do.call(rbind, .)


estimate_comparison_df <- lapply(c("alpha", "rho"), function(par) {
  
  estimate_comparison_df %>% 
    filter(parameter == par) %>% 
    melt() %>% 
    select(variable, value) %>% 
    set_colnames(c("variable", par))
  
}) %>% do.call(cbind, .) %>% 
  set_colnames(c("dummy", "alpha", "variable", "rho")) %>% 
  select(variable, rho, alpha) %>%
  mutate(index = rep(1:100, 3))


png("figures/manual_computations/compare_estimates.png")
p <- estimate_comparison_df %>% 
  ggplot() +
  geom_point(aes(x = rho, y = alpha, color = variable)) +
  geom_line(aes(x = rho, y = alpha, group = index)) +
  scale_color_manual(labels = c("Stan", "Truth", "Manual MAP"),
                     values = c("red", "black", "blue"))
print(p)
dev.off()




## DUMP ***** ####
dat <- single_simulated_data[[rho_true]][[alpha_true]] %>% 
  filter((x%%1==0))
dat <- list(N = nrow(dat), y = dat$y, x = dat$x)




# likelihood

likelihood_val <- data.frame(rho = rep(rho_lh_grid, length(alpha_lh_grid)),
                             alpha = rep(alpha_lh_grid, each = length(rho_lh_grid)),
                             value = NA)

for(alpha in alpha_lh_grid) {
  
  for(rho in rho_lh_grid) {
    
    condition <- rho == likelihood_val$rho & alpha == likelihood_val$alpha
    
    
    covM <- matrix(NA, length(dat$x), length(dat$x))
    
    for(i in 1:length(dat$x)) {
      for(j in 1:length(dat$x)) {
        
        covM[i, j] <- alpha^2*exp(-abs(i-j)/(rho^2)) + ifelse(i == j, measurement_error^2, 0)
        
      }
    }
    
    
    # likelihood_val[condition, "value"] <- dmvnorm(x = dat$y, sigma = covM)
    
  }
  
}


posterior_val <- cbind(likelihood_val, prior = prior$value) %>% 
  set_colnames(c("rho", "alpha", "likelihood", "prior")) %>% 
  mutate(posterior = likelihood*prior) %>% 
  mutate(posterior = posterior/sum(posterior)) %>%
  mutate(log_likelihood = log(likelihood)) %>% 
  melt(id.vars = c("alpha", "rho"))

# ## PC Prior ******* 
# 
# alpha0 <- 1
# rho0 <- 10
# tail <- .5
# 
# lambda1 = -log(tail)*sqrt(2*rho0)
# lambda2 = -log(tail)/alpha0
# 
# pc_prior_val <- data.frame(rho = rep(1:100/10, 50),
#                            alpha = rep(1:50/10, each = 100),
#                            value = NA)
# 
# for(alpha in 1:50/10) {
#   
#   for(rho in 1:100/10) {
#     
#     condition <- rho == pc_prior_val$rho & alpha == pc_prior_val$alpha
#     
#     pc_prior_val[condition, "value"] <- lambda1*lambda2*(2*rho)^(-3/2)*exp(-lambda1*sqrt(2*rho) - lambda2*alpha)
#     
#     
#   }
#   
# }
# 
# 
# 
# pc_prior_val %>% 
#   ggplot(aes(x = rho, y = alpha,  size = value)) +
#   geom_point()
# 
# 
