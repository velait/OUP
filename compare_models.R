# Compare different models




## Data **************** ####

n_series <- 5
intervals <- 1:30

# # lambda ~ logNormal
# lambda_mean <- -.5
# lambda_sd <- 0.1
# lambda <- restricted_rlnorm(n_series, meanlog = lambda_mean, sdlog = lambda_sd)

# lambda ~ inv_gamma
lambda_mean <- 20
lambda_sd <- 10
lambda <- oup_invG_lambda(n_series, shape = lambda_mean, scale = lambda_sd)


# sigma ~ logNormal
sigma_mean <- .2
sigma_sd <- .1
sigma <- rnorm(n_series, mean = sigma_mean, sd = sigma_sd)


# mu ~ Normal
mu_mean <- 0
mu_sd <- .25
mu <- rnorm(n_series, mean = mu_mean, sd = mu_sd)

# Generate set of series
diff_compare_data <- generate_student_set(n_series = n_series,
                                          student_df = 7,
                                          mu = mu,
                                          sigma =  sigma,
                                          lambda = lambda,
                                          intervals = intervals,
                                          seed = 11235)

hyper_parameter_values <- c(mu_mean = mu_mean,
                            mu_sd = mu_sd, 
                            sigma_mean = sigma_mean,
                            sigma_sd = sigma_sd,
                            lambda_mean = lambda_mean,
                            lambda_sd = lambda_sd)

## Models ************** ####
pooled_student_t_oup <- stan_model("stan_models/pooled_student_t_oup.stan")
non_pooled_student_t_oup <- stan_model("stan_models/non_pooled_student_t_oup.stan")
hierarchical_student_t_oup <- stan_model("stan_models/hierarchical_student_t_oup.stan")


## Samples ************* ####

pooled_samples <- sampling(pooled_student_t_oup,
                                diff_compare_data,
                                iter=iter,
                                chains=chains,
                                init=1)

non_pooled_samples <- sampling(non_pooled_student_t_oup,
                                    diff_compare_data,
                                    iter=iter,
                                    chains=chains, 
                                    init=1)

partially_pooled_samples <- sampling(hierarchical_student_t_oup,
                                          diff_compare_data,
                                          iter=iter,
                                          chains=chains,
                                          init=1)

save(pooled_samples, 
     non_pooled_samples,
     partially_pooled_samples,
     lambda, 
     mu,
     sigma,
     intervals, 
     n_series,
     file = "results/samples_5_30.Rds")

## Results ************* ####


# Stan samples to a list
compare_models_list <- list(pooled_samples = pooled_samples,
                                  non_pooled_samples = non_pooled_samples, 
                                  partially_pooled_samples = partially_pooled_samples)






compare_summary <- lapply(names(compare_models_list), function(x) {
   summary(compare_models_list[[x]])$summary[grep("lambda\\[|sigma\\[|mu\\[", rownames(summary(compare_models_list[[x]])$summary)), c("25%", "50%", "75%")] %>%
    as.data.frame() %>%
    rownames_to_column(var="parameter") %>%
    set_colnames(c("parameter", "lower", "mode", "upper"))
  
}) %>% set_names(names(compare_models_list))


# Add pooled results
compare_summary[["pooled_samples"]] <- summary(compare_models_list[["pooled_samples"]])$summary[grep("lambda|sigma|mu|lp__", rownames(summary(compare_models_list[["pooled_samples"]])$summary)), c("25%", "50%", "75%")] %>%
  as.data.frame() %>%
  rownames_to_column(var="parameter") %>%
  set_colnames(c("parameter", "lower", "mode", "upper"))


# Hyper parameters
hyper_parameter_summary <-  summary(partially_pooled_samples)$summary[hyper_parameter_values %>%
                                                                        names, c("25%", "50%", "75%")]



## Compare learned prior to simulation priors
hyper_prior_plots <- lapply(c("lambda", "mu", "sigma"), function(par) {
  
  if(par == "mu") {
    
    prior <- dnorm(seq(from=-2, to=2, length.out = 100),
                   mu_mean,
                   mu_sd)
    
    posterior <- dnorm(seq(from=-2, to=2, length.out = 100),
                       hyper_parameter_summary["mu_mean", "50%"],
                       hyper_parameter_summary["mu_sd", "50%"])
    
    both <- cbind(Prior=prior, Posterior=posterior) %>%
      melt %>%
      mutate(x = rep(seq(from=-4, to=4, length.out = 100),2))
    
    p <-  ggplot(both, aes(x=x, y=value, color=Var2)) +
      geom_line() +
      guides(linetype=FALSE) +
      labs(y="", title = par) +
      xlab(expression(mu)) +
      guides(color=guide_legend(title="Model"))
    
    
    
    
  } else if(par == "lambda") {
    
    prior <- dlnorm(seq(from=0, to=1, length.out = 100),
                   lambda_mean,
                   lambda_sd)
    
    posterior <- dlnorm(seq(from=0, to=1, length.out = 100),
                       hyper_parameter_summary["lambda_mean", "50%"],
                       hyper_parameter_summary["lambda_sd", "50%"])
    
    both <- cbind(Prior=prior, Posterior=posterior) %>%
      melt %>%
      mutate(x = rep(seq(from=0, to=1, length.out = 100),2))
    
    p <-  ggplot(both, aes(x=x, y=value, color=Var2)) +
      geom_line() +
      guides(linetype=FALSE) +
      labs(y="", title = par) +
      xlab(expression(lambda)) +
      guides(color=guide_legend(title="Model"))
    
  } else if(par == "sigma") {
    
    prior <- dlnorm(seq(from=0, to=1, length.out = 100),
                    sigma_mean,
                    sigma_sd)
    
    posterior <- dlnorm(seq(from=0, to=1, length.out = 100),
                        hyper_parameter_summary["sigma_mean", "50%"],
                        hyper_parameter_summary["sigma_sd", "50%"])
    
    both <- cbind(Prior=prior, Posterior=posterior) %>%
      melt %>%
      mutate(x = rep(seq(from=0, to=1, length.out = 100),2))
    
    p <-  ggplot(both, aes(x=x, y=value, color=Var2)) +
      geom_line() +
      guides(linetype=FALSE) +
      labs(y="", title = par) +
      xlab(expression(sigma)) +
      guides(color=guide_legend(title="Model"))
    
    
  } 

  
  
  
}) 


## Plots *************** ####
diff_compare_plots <- list()
diff_partially_pooled_errorbar_plots <- list()
diff_no_pooling_errorbar_plots <- list()
diff_cross_plots <- list()

for(par in c("lambda", "sigma", "mu")) {
  
  pooled_estimate <- (compare_summary[["pooled_samples"]] %>% 
                        filter(parameter==par) %>%
                        mutate(model="Complete pooling", series=NA))[, "mode"] 
  
  # pooled_estimate[1:40, ] <- compare_summary[["pooled_samples"]] %>% filter(parameter==par) %>% mutate(model="Complete pooling", series=NA)
  # 
  # pooled_estimate$series <- 1:40
  
  ord <- diff_compare_data[[paste0(par, "_values")]] %>% order
  
  
  
  partial_estimate <- compare_summary[["partially_pooled_samples"]] %>%
    filter(str_detect(parameter, par)) %>%
    mutate(series = str_extract(parameter, "\\d+") %>%
             as.numeric(), model="Partially pooled")
  
  partial_estimate <- partial_estimate[match(ord, partial_estimate$series),] %>%
    mutate(ord = 1:n_series)
  
  non_pooled_estimate <- compare_summary[["non_pooled_samples"]] %>% 
    filter(str_detect(parameter, par)) %>%
    mutate(series = str_extract(parameter, "\\d+") %>%
             as.numeric(), model="Non pooled")
  
  non_pooled_estimate <- non_pooled_estimate[match(ord, non_pooled_estimate$series),] %>% mutate(ord = 1:n_series)
  
  estimates <- rbind(partial_estimate, non_pooled_estimate) %>% dplyr::select(-parameter)
  
  p <- ggplot() +
    geom_point(data=estimates, aes(x=ord, y=mode, shape=model))  +
    labs(x="Time series ID", y="Posterior estimate") + 
    guides(shape=guide_legend("")) +
    geom_hline(yintercept=pooled_estimate, linetype = "dashed") +
    geom_line(data=diff_compare_data[[paste0(par, "_values")]]%>% sort() %>% as_tibble(), aes(y=value, x=1:n_series)) +
    scale_shape_manual(values=c(1, 16), labels=c("No pooling", "Partial pooling")) +
    theme(legend.position="top")
  
  diff_compare_plots[[par]] <- p
  
  
  q <- ggplot()  + 
    labs(x="Series", y="Posterior estimate", title=par) +
    guides(shape=guide_legend("")) +
    geom_line(data=diff_compare_data[[paste0(par, "_values")]]%>% sort() %>% as_tibble(), aes(y=value, x=1:n_series)) +
    scale_shape_manual(values=c(1, 16)) + 
    geom_errorbar(data=estimates %>% filter(model=="Partially pooled"), aes(x=ord, ymin=lower, ymax=upper)) 
  
  diff_partially_pooled_errorbar_plots[[par]] <- q
  
  pp <- ggplot()  + 
    labs(x="Series", y="Posterior estimate", title=par) +
    guides(shape=guide_legend("")) +
    geom_line(data=diff_compare_data[[paste0(par, "_values")]]%>% sort() %>% as_tibble(), aes(y=value, x=1:n_series)) +
    scale_shape_manual(values=c(1, 16)) + 
    geom_errorbar(data = estimates, aes(x=ord, ymin=lower, ymax=upper, color = model)) 
  
  
  diff_no_pooling_errorbar_plots[[par]] <- pp
  
  # Cross plots
  
  m <- max(c(diff_compare_data[[paste0(par, "_values")]]%>% sort(), partial_estimate$mode))
  
  mi <- min(c(diff_compare_data[[paste0(par, "_values")]]%>% sort(), partial_estimate$mode))
  vals <- diff_compare_data[[paste0(par, "_values")]]%>% sort()
  corr <- cor(vals, partial_estimate$mode, method = "spearman")
  cross_p <- ggplot(cbind(vals, partial_estimate), aes(x=mode, y=vals)) + geom_point() + geom_abline(intercept = 0, slope=1) + scale_x_continuous(limits = c(mi, m))+ scale_y_continuous(limits = c(mi, m)) + labs(title=paste(par, "correlation: ", round(corr, 2)))
  
  diff_cross_plots[[par]] <- cross_p
  
}

