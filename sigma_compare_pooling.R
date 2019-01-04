#### Pooled vs. partially pooled vs. non pooled models
#### Samples sizes comparable to those in Ravel et al used.

#### DATA ####

mu_mu <- 0
mu_sigma <- 0.001
sigma_mu <- .5
sigma_sigma <- 0.001
lambda_alpha <- 20
lambda_beta <- 10

compare_n_series <- 15
mu_val <- rnorm(compare_n_series, mu_mu, mu_sigma)
sigma_val <- rnorm(compare_n_series, sigma_mu, sigma_sigma)
lambda_val <- oup_invG_lambda(compare_n_series, shape=lambda_alpha, scale=lambda_beta)
kappa_val <- (sigma_val^2)/(2*lambda_val)

hyper_parameter_values <- c(mu_mu=mu_mu,
                            mu_sigma=mu_sigma,
                            sigma_mu=sigma_mu,
                            sigma_sigma=sigma_sigma,
                            lambda_alpha=lambda_alpha,
                            lambda_beta=lambda_beta,
                            inv_gamma_mode=lambda_beta/(lambda_alpha+1),
                            inv_gamma_variance=lambda_beta^2/((lambda_alpha-2)*(lambda_alpha-1)^2))

diff_compare_data <- generate_student_set(n_series = compare_n_series,
                                          student_df = 7,
                                          mu = mu_val,
                                          sigma =  sigma_val,
                                          kappa = kappa_val,
                                          lambda = lambda_val,
                                          intervals = 1:20,
                                          seed = 11235)

# diff_compare_data <- generate_student_set(n_series = compare_n_series, student_df = 7, mu = mu_val, sigma =  sigma_val, lambda = lambda_val, intervals = 1:10)

diff_compare_data[["lambda_values"]] %>% hist(breaks=10)



#### MODELS ####

# pooled_student_t_oup <- stan_model("pooled_student_t_oup.stan")
non_pooled_student_t_oup <- stan_model("stan_models/kappa_hierarchical_student_t_oup.stan")
hierarchical_student_t_oup <- stan_model("stan_models/hierarchical_student_t_oup.stan")


#### SAMPLES ####
# diff_pooled_samples <- sampling(pooled_student_t_oup, diff_compare_data, iter=iter, chains=chains, init=1)

diff_non_pooled_samples <- sampling(non_pooled_student_t_oup, diff_compare_data, iter=iter, chains=chains, init=1)

diff_partially_pooled_samples <- sampling(hierarchical_student_t_oup, diff_compare_data, iter=iter, chains=chains, init=1)



#### RESULTS ####

# Stan samples to a list
diff_compare_samples_list <- list(pooled_samples=diff_pooled_samples, non_pooled_samples=diff_non_pooled_samples, partially_pooled_samples=diff_partially_pooled_samples)

# save(diff_compare_samples_list, file="ravel_simulation_samples")

# Parameters
diff_compare_summary <- list()
for(x in names(diff_compare_samples_list)) {
  diff_compare_summary[[x]] <- summary(diff_compare_samples_list[[x]])$summary[grep("lambda\\[|sigma\\[|mu\\[", rownames(summary(diff_compare_samples_list[[x]])$summary)), c("25%", "50%", "75%")] %>% as.data.frame() %>% rownames_to_column(var="parameter") %>% set_colnames(c("parameter", "lower", "mode", "upper"))
}

diff_compare_summary[["pooled_samples"]] <- summary(diff_compare_samples_list[["pooled_samples"]])$summary[grep("lambda|sigma|mu|lp__", rownames(summary(diff_compare_samples_list[["pooled_samples"]])$summary)), c("25%", "50%", "75%")] %>% as.data.frame() %>% rownames_to_column(var="parameter") %>% set_colnames(c("parameter", "lower", "mode", "upper"))

# Hyper paramters
diff_hyper_parameter_summary <-  summary(diff_partially_pooled_samples)$summary[hyper_parameter_values %>% names, c("25%", "50%", "75%")]

## Compare learned prior to simulation priors

ravel_hyper_prior_plots <- list()
for(par in c("lambda", "mu", "sigma")) {
  
  if(grepl("mu", par)) {
    par <- "mu"
    prior <- dnorm(seq(from=-4, to=4, length.out = 100), 0, 1)
    posterior <- dnorm(seq(from=-4, to=4, length.out = 100), diff_hyper_parameter_summary["mu_mu", "50%"], diff_hyper_parameter_summary["mu_sigma", "50%"])
    
    both <- cbind(Prior=prior, Posterior=posterior) %>% melt %>% mutate(x = rep(seq(from=-4, to=4, length.out = 100),2))
    
    ravel_hyper_prior_plots[[par]] <-  ggplot(both, aes(x=x, y=value, linetype=Var2)) + geom_line() + guides(linetype=FALSE) +labs(y="") + xlab(expression(mu))
    
    
    ggplot()
    
  } else if(grepl("lambda", par)) {
    par <- 'lambda'
    prior <- dinvgamma(seq(from=0, to=1, length.out = 100), 6, 4)
    posterior <- dinvgamma(seq(from=0, to=1, length.out = 100), diff_hyper_parameter_summary["lambda_alpha", "50%"], diff_hyper_parameter_summary["lambda_beta", "50%"])
    
    both <- cbind(Prior=prior, Posterior=posterior) %>% melt %>% mutate(x = rep(seq(from=0, to=1, length.out = 100),2))
    
    ravel_hyper_prior_plots[[par]] <-  ggplot(both, aes(x=x, y=value, linetype=Var2)) + geom_line() + guides(linetype=FALSE) +labs(y="Density") + xlab(expression(lambda))
    
  } else if(grepl("sigma", par)) {
    
    par <- 'sigma'
    prior <- dnorm(seq(from=0, to=10, length.out = 100), 3, 1)
    posterior <- dnorm(seq(from=0, to=10, length.out = 100), diff_hyper_parameter_summary["sigma_mu", "50%"], diff_hyper_parameter_summary["sigma_sigma", "50%"])
    
    both <- cbind(Prior=prior, Posterior=posterior) %>% melt %>% mutate(x = rep(seq(from=0, to=10, length.out = 100),2))
    
    ravel_hyper_prior_plots[[par]] <-  ggplot(both, aes(x=x, y=value, linetype=Var2)) + geom_line() + guides(linetype=FALSE) +labs(y="") + xlab(expression(sigma))
    
  }
  
  # dat <- tibble(prior=prior, posterior=posterior) %>% melt
  # 
  # ravel_hyper_prior_plots[[par]] <-  ggplot(dat) +
  #   geom_density(aes(value, fill=variable), alpha=0.1) 
  
  
  
  
}


#### PLOT ALL ####

diff_compare_plots <- list()
diff_partially_pooled_errorbar_plots <- list()
diff_no_pooling_errorbar_plots <- list()
diff_cross_plots <- list()

for(par in c("lambda", "sigma", "mu")) {
  
  pooled_estimate <- (diff_compare_summary[["pooled_samples"]] %>% filter(parameter==par) %>% mutate(model="Complete pooling", series=NA))[, "mode"] 
  
  # pooled_estimate[1:40, ] <- compare_summary[["pooled_samples"]] %>% filter(parameter==par) %>% mutate(model="Complete pooling", series=NA)
  # 
  # pooled_estimate$series <- 1:40
  
  ord <- diff_compare_data[[paste0(par, "_values")]] %>% order
  
  # simulation_value <- oup_simulation_parameters[par]
  
  partial_estimate <- diff_compare_summary[["partially_pooled_samples"]] %>% filter(str_detect(parameter, par)) %>% mutate(series = str_extract(parameter, "\\d+") %>% as.numeric(), model="Partially pooled")
  
  partial_estimate <- partial_estimate[match(ord, partial_estimate$series),] %>% mutate(ord = 1:compare_n_series)
  
  non_pooled_estimate <- diff_compare_summary[["non_pooled_samples"]] %>% filter(str_detect(parameter, par)) %>% mutate(series = str_extract(parameter, "\\d+") %>% as.numeric(), model="No pooling")
  
  non_pooled_estimate <- non_pooled_estimate[match(ord, non_pooled_estimate$series),] %>% mutate(ord = 1:compare_n_series)
  
  estimates <- rbind(partial_estimate, non_pooled_estimate) %>% dplyr::select(-parameter)
  
  p <- ggplot() +
    geom_point(data=estimates, aes(x=ord, y=mode, shape=model))  +
    labs(x="Time series ID", y="Posterior estimate") + 
    guides(shape=guide_legend("")) +
    geom_hline(yintercept=pooled_estimate, linetype = "dashed") +
    geom_line(data=diff_compare_data[[paste0(par, "_values")]]%>% sort() %>% as_tibble(), aes(y=value, x=1:compare_n_series)) +
    scale_shape_manual(values=c(1, 16), labels=c("No pooling", "Partial pooling")) +
    theme(legend.position="top")
  
  diff_compare_plots[[par]] <- p
  
  
  q <- ggplot()  + 
    labs(x="Series", y="Posterior estimate", title=par) +
    guides(shape=guide_legend("")) +
    geom_line(data=diff_compare_data[[paste0(par, "_values")]]%>% sort() %>% as_tibble(), aes(y=value, x=1:compare_n_series)) +
    scale_shape_manual(values=c(1, 16)) + 
    geom_errorbar(data=estimates %>% filter(model=="Partially pooled"), aes(x=ord, ymin=lower, ymax=upper)) 
  
  diff_partially_pooled_errorbar_plots[[par]] <- q
  
  qq <- ggplot()  + 
    labs(x="Series", y="Posterior estimate", title=par) +
    guides(shape=guide_legend("")) +
    geom_line(data=diff_compare_data[[paste0(par, "_values")]]%>% sort() %>% as_tibble(), aes(y=value, x=1:compare_n_series)) +
    scale_shape_manual(values=c(1, 16)) + 
    geom_errorbar(data=estimates %>% filter(model=="No pooling"), aes(x=ord, ymin=lower, ymax=upper)) 
  
  diff_no_pooling_errorbar_plots[[par]] <- qq
  
  
  # Cross plots
  
  m <- max(c(diff_compare_data[[paste0(par, "_values")]]%>% sort(), partial_estimate$mode))
  
  mi <- min(c(diff_compare_data[[paste0(par, "_values")]]%>% sort(), partial_estimate$mode))
  vals <- diff_compare_data[[paste0(par, "_values")]]%>% sort()
  corr <- cor(vals, partial_estimate$mode, method = "spearman")
  cross_p <- ggplot(cbind(vals, partial_estimate), aes(x=mode, y=vals)) + geom_point() + geom_abline(intercept = 0, slope=1) + scale_x_continuous(limits = c(mi, m))+ scale_y_continuous(limits = c(mi, m)) + labs(title=paste(par, "correlation: ", round(corr, 2)))
  
  diff_cross_plots[[par]] <- cross_p
  
}




#### PLOT IQR ####

diff_error_plots <- list()
diff_IQR_plots <- list()
for(par in c("lambda", "sigma", "mu")) {
  
  real_values <- diff_compare_data[[paste0(par, "_values")]]
  ord <- real_values  %>% order
  
  pooled_estimate <- (diff_compare_summary[["pooled_samples"]] %>% filter(parameter==par) %>% mutate(model="Complete pooling", series=NA))
  
  pooled_estimate[1:compare_n_series, ] <- diff_compare_summary[["pooled_samples"]] %>% filter(parameter==par) %>% mutate(model="Complete pooling", series=NA) 
  
  pooled_estimate <- pooled_estimate %>% mutate(real_value=real_values, series=1:compare_n_series)
  
  # pooled_estimate$series <- 1:80
  
  pooled_estimate <- pooled_estimate[match(ord, pooled_estimate$series),] %>% mutate(ord = 1:compare_n_series)
  
  
  
  # simulation_value <- oup_simulation_parameters[par]
  
  
  # PARTIAL
  partial_estimate <- diff_compare_summary[["partially_pooled_samples"]] %>% filter(str_detect(parameter, par)) %>% mutate(series = str_extract(parameter, "\\d+") %>% as.numeric(), model="Partially pooled", real_value=real_values)
  
  partial_estimate <- partial_estimate[match(ord, partial_estimate$series),] %>% mutate(ord = 1:compare_n_series)
  
  # NON
  non_pooled_estimate <- diff_compare_summary[["non_pooled_samples"]] %>% filter(str_detect(parameter, par)) %>% mutate(series = str_extract(parameter, "\\d+") %>% as.numeric(), model="No pooling", real_value=real_values)
  
  non_pooled_estimate <- non_pooled_estimate[match(ord, non_pooled_estimate$series),] %>% mutate(ord = 1:compare_n_series)
  
  # estimates <- rbind(pooled_estimate, partial_estimate, non_pooled_estimate) %>% dplyr::select(-parameter) %>% mutate(error=mode-real_value, IQR=upper-lower) 
  
  estimates <- rbind(pooled_estimate, partial_estimate, non_pooled_estimate) %>% dplyr::select(-parameter) %>% mutate(error=abs(mode-real_value)/real_value, IQR=upper-lower) 
  
  estimates$model <-  factor(estimates$model, levels = c("No pooling", "Partially pooled", "Complete pooling"), ordered = TRUE)
  
  p <- ggplot(estimates)+ geom_hline(yintercept = 0, linetype="dashed") +
    geom_boxplot(aes(x=model, y=error), width=0.5) +
    labs(x="", y="Error") + scale_x_discrete(labels=c("No pooling" = "None", "Partially pooled" = "Partial", "Complete pooling" = "Complete")) 
  
  diff_error_plots[[par]] <- p
  
  q <- ggplot(estimates) + 
    geom_boxplot(aes(x=model, y=IQR), width=0.5) +
    labs(x="",  y="50% IQR width")+ scale_x_discrete(labels=c("No pooling" = "None", "Partially pooled" = "Partial", "Complete pooling" = "Complete"))
  
  diff_IQR_plots[[par]] <- q
}


#### COWPLOT ####
bottom_row <- plot_grid(diff_error_plots[['lambda']], diff_IQR_plots[['lambda']], labels=c("B", "C"))

plot_grid(diff_compare_plots[['lambda']], bottom_row, labels=c("A", ""), ncol = 1, rel_heights = c(1, 1))



