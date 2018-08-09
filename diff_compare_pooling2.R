#### Pooled vs. partially pooled vs. non pooled models

#### DATA ####

oup_invG_lambda <- function(n=compare_n_series, shape=2, scale=0.5) {
  vec <- c()
  for(i in 1:n) {
    
    x <- 1
    while(x >= 1) {
      x <- MCMCpack::rinvgamma(1, shape = 2, scale = 0.5)
    }
    
    vec[i] <- x
  }
  
  return(vec)
}

oup_G_lambda <- function(n=compare_n_series, shape=2, scale=4) {
  vec <- c()
  for(i in 1:n) {
    
    x <- 1
    while(x >= 1) {
      x <- rgamma(1, shape = shape, scale = scale)
    }
    
    vec[i] <- x
  }
  
  return(vec)
}

oup_normal_lambda <- function(n=compare_n_series, shape=0.4, scale=0.2) {
  vec <- c()
  for(i in 1:n) {
    
    x <- 1
    while(x >= 1 | x <= 0) {
      x <- rnorm(1, shape, scale)
    }
    
    vec[i] <- x
  }
  
  return(vec)
}

compare_n_series <- 40
# diff_compare_data <- generate_student_set(n_series = compare_n_series, student_df = 7, mu = rnorm(compare_n_series, mu_mu, mu_sigma), sigma = rnorm(compare_n_series, sigma_mu, sigma_sigma), lambda = oup_invG_lambda(compare_n_series, shape=2, scale=0.5), intervals = 1:20)

mu_val <- rep(5, compare_n_series)
sigma_val <- rep(0.2, compare_n_series)
lambda_val <- oup_invG_lambda(compare_n_series, 4, 0.4)


diff_compare_data <- generate_student_set(n_series = compare_n_series, student_df = 7, mu = mu_val, sigma =  sigma_val, lambda = lambda_val, intervals = 1:10)



#### MODELS ####
pooled_student_t_oup <- stan_model("pooled_student_t_oup.stan")
non_pooled_student_t_oup <- stan_model("non_pooled_student_t_oup.stan")
hierarchical_student_t_oup <- stan_model("hierarchical_student_t_oup.stan")


#### SAMPLES ####
diff_pooled_samples <- sampling(pooled_student_t_oup, diff_compare_data, iter=iter, chains=chains, init=1)

diff_non_pooled_samples <- sampling(non_pooled_student_t_oup, diff_compare_data, iter=iter, chains=chains, init=1)

diff_partially_pooled_samples <- sampling(hierarchical_student_t_oup, diff_compare_data, iter=iter, chains=chains, init=0.5)



#### RESULTS ####

# get hyper parameter summaries in list
diff_compare_samples_list <- list(pooled_samples=diff_pooled_samples, non_pooled_samples=diff_non_pooled_samples, partially_pooled_samples=diff_partially_pooled_samples) 

# save(diff_compare_samples_list, file="diff_compare_samles_list")

diff_compare_summary <- list()
for(x in names(diff_compare_samples_list)) {
  diff_compare_summary[[x]] <- summary(diff_compare_samples_list[[x]])$summary[grep("lambda\\[|sigma\\[|mu\\[", rownames(summary(diff_compare_samples_list[[x]])$summary)), c("25%", "50%", "75%")] %>% as.data.frame() %>% rownames_to_column(var="parameter") %>% set_colnames(c("parameter", "lower", "mode", "upper"))
}

diff_compare_summary[["pooled_samples"]] <- summary(diff_compare_samples_list[["pooled_samples"]])$summary[grep("lambda|sigma|mu|lp__", rownames(summary(diff_compare_samples_list[["pooled_samples"]])$summary)), c("25%", "50%", "75%")] %>% as.data.frame() %>% rownames_to_column(var="parameter") %>% set_colnames(c("parameter", "lower", "mode", "upper"))



#### PLOT ALL ####

diff_compare_plots <- list()
diff_partially_pooled_errorbar_plots <- list()
diff_no_pooling_errorbar_plots <- list()

for(par in c("lambda", "sigma", "mu")) {
  
  pooled_estimate <- (diff_compare_summary[["pooled_samples"]] %>% filter(parameter==par) %>% mutate(model="Complete pooling", series=NA))[, "mode"] 
  
  # pooled_estimate[1:40, ] <- compare_summary[["pooled_samples"]] %>% filter(parameter==par) %>% mutate(model="Complete pooling", series=NA)
  # 
  # pooled_estimate$series <- 1:40
  
  ord <- diff_compare_data[[paste0(par, "_values")]] %>% order
  
  simulation_value <- oup_simulation_parameters[par]
  
  partial_estimate <- diff_compare_summary[["partially_pooled_samples"]] %>% filter(str_detect(parameter, par)) %>% mutate(series = str_extract(parameter, "\\d+") %>% as.numeric(), model="Partially pooled")
  
  partial_estimate <- partial_estimate[match(ord, partial_estimate$series),] %>% mutate(ord = 1:compare_n_series)
  
  non_pooled_estimate <- diff_compare_summary[["non_pooled_samples"]] %>% filter(str_detect(parameter, par)) %>% mutate(series = str_extract(parameter, "\\d+") %>% as.numeric(), model="No pooling")
  
  non_pooled_estimate <- non_pooled_estimate[match(ord, non_pooled_estimate$series),] %>% mutate(ord = 1:compare_n_series)
  
  estimates <- rbind(partial_estimate, non_pooled_estimate) %>% dplyr::select(-parameter)
  
  p <- ggplot() +
    geom_point(data=estimates, aes(x=ord, y=mode, shape=model))  +
    labs(x="Series", y="Posterior estimate", title=paste0("\\", par)) + 
    guides(shape=guide_legend("")) +
    geom_hline(yintercept=pooled_estimate, linetype = "dashed") +
    geom_line(data=diff_compare_data[[paste0(par, "_values")]]%>% sort() %>% as_tibble(), aes(y=value, x=1:compare_n_series)) +
    scale_shape_manual(values=c(1, 16))
  
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
  
  estimates <- rbind(pooled_estimate, partial_estimate, non_pooled_estimate) %>% dplyr::select(-parameter) %>% mutate(error=mode-real_value, IQR=upper-lower)
  
  p <- ggplot(estimates) +
    geom_boxplot(aes(x=model, y=error)) +
    labs(y="Posterior mode - simulation value", x="Model", title=par)
  
  diff_error_plots[[par]] <- p
  
  q <- ggplot(estimates) + 
    geom_boxplot(aes(x=model, y=IQR)) +
    labs(x="Model", y="50% interquartile range", title=par)
  
  diff_IQR_plots[[par]] <- q
}


#### SAVE PLOTS ####



  # setEPS()
  # postscript(paste0("plots/diff_compare_plots_lambda.eps"))
  # diff_compare_plots[["lambda"]]
  # dev.off()
  # 
  # setEPS()
  # postscript(paste0("plots/diff_compare_plots_mu.eps"))
  # diff_compare_plots[["mu"]]
  # dev.off()
  # 
  # setEPS()
  # postscript(paste0("plots/diff_compare_plots_sigma.eps"))
  # diff_compare_plots[["sigma"]]
  # dev.off()
  # 
  # 
  # 
  # setEPS()
  # postscript(paste0("plots/diff_partially_pooled_errorbar_plots_lambda.eps"))
  # diff_partially_pooled_errorbar_plots[["lambda"]]
  # dev.off()
  # 
  # setEPS()
  # postscript(paste0("plots/diff_partially_pooled_errorbar_plots_plots_mu.eps"))
  # diff_partially_pooled_errorbar_plots[["mu"]]
  # dev.off()
  # 
  # setEPS()
  # postscript(paste0("plots/diff_partially_pooled_errorbar_plots_sigma.eps"))
  # diff_partially_pooled_errorbar_plots[["sigma"]]
  # dev.off()
  # 
  # 
  # setEPS()
  # postscript(paste0("plots/diff_no_pooling_errorbar_plots_lambda.eps"))
  # diff_no_pooling_errorbar_plots[["lambda"]]
  # dev.off()
  # 
  # setEPS()
  # postscript(paste0("plots/diff_no_pooling_errorbar_plots_mu.eps"))
  # diff_no_pooling_errorbar_plots[["mu"]]
  # dev.off()
  # 
  # setEPS()
  # postscript(paste0("plots/diff_no_pooling_errorbar_plots_sigma.eps"))
  # diff_no_pooling_errorbar_plots[["sigma"]]
  # dev.off()
  # 
  # 
  # setEPS()
  # postscript(paste0("plots/diff_error_plots_lambda.eps"))
  # diff_error_plots[["lambda"]]
  # dev.off()
  # 
  # setEPS()
  # postscript(paste0("plots/diff_error_plots_mu.eps"))
  # diff_error_plots[["mu"]]
  # dev.off()
  # 
  # setEPS()
  # postscript(paste0("plots/diff_error_plots_sigma.eps"))
  # diff_error_plots[["sigma"]]
  # dev.off()
  # 
  # 
  # 
  # 
  # setEPS()
  # postscript(paste0("plots/diff_IQR_plots_lambda.eps"))
  # diff_IQR_plots[["lambda"]]
  # dev.off()
  # 
  # setEPS()
  # postscript(paste0("plots/diff_IQR_plots_mu.eps")) 
  # diff_IQR_plots[["mu"]]
  # dev.off()
  # 
  # setEPS()
  # postscript(paste0("plots/diff_IQR_plots_sigma.eps"))
  # diff_IQR_plots[["sigma"]]
  # dev.off()



