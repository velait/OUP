library(tidyverse)
library(microbiome)
library(magrittr)
library(reshape2)
## Functions ****************************** ####

# Function to get predicted values for a given taxon
bimodality_premonition <- function(delta_df, taxon = "Dialister") {
  
  change_points <- c("Anaerofustis" = 0.4, 
                     "Aneurinibacillus" = 0.3, 
                     "Aquabacterium" = 0.5, 
                     "Bacteroides intestinalis et rel." = .5, 
                     "Bacteroides fragilis et rel." = 3.1,
                     "Burkholderia" = 0.5, 
                     "Dialister" = 3.25, 
                     "Leminorella" = 0.5, 
                     "Prevotella melaninogenica et rel." = 4.5, 
                     "Prevotella oralis et rel." = 4, 
                     "Serratia" = 0.5, 
                     "Uncultured Bacteroidetes" = 0.5, 
                     "Uncultured Clostridiales II" = 3.4, 
                     "Uncultured Selenomonadaceae" = 0.5, 
                     "Wissella et rel." = 0.45)
  
  change_point <- change_points[taxon]
  
  # otu_abundance <- data.frame(abundance = t(abundances(atlas1006_log))[, taxon])
  
  otu_df <- delta_df %>% 
    filter(otu == taxon) %>% 
    mutate(end = start + delta_abundance) %>% 
    mutate(start_mode = ifelse(start < change_point, "low", "high")) %>% 
    mutate(end_mode = ifelse(start + delta_abundance < change_point, "low", "high")) %>% 
    mutate(change = ifelse(start_mode != end_mode, "yes", "no"),
           change_to_high = ifelse(start_mode == "low" & end_mode == "high", "yes", "no"), 
           change_to_low = ifelse(start_mode == "high" & end_mode == "low", "yes", "no"))
  
  
  # Fit polynomial to rates
  smooth_rate <- loess(unit_rate ~ start, otu_df)
  
  
  # Get predicted rate of change for each observation which isn't the subject's last observation
  # x_seq <- seq(from = range(otu_df$start)[1], to = range(otu_df$start)[2], by = 0.01)
  # 
  # rate_prediction <- predict(smooth_rate, newdata =  data.frame(start = x_seq), se = TRUE)
  # rate_prediction_df <- data.frame(rate = rate_prediction$fit, x = x_seq)
  
  # get predicted end point
  # otu_df$end_prediction <- otu_df$start + rate_prediction$fit*otu_df$delta_t # SIMPLE
  
  otu_df$end_prediction <- lapply(1:nrow(otu_df), function(i) {
    
    my_rate <- otu_df[i, c("delta_t", "start")]
    
    # Number of time steps
    stepsize <- 0.1
    n_dt <- my_rate$delta_t/stepsize
    
    # x <- rep(NA, n_dt)
    # x[1] <- my_rate$start
    x <- my_rate$start
    
    for(t in 1:n_dt) {
      # x[t+1] <- x[t] + stepsize*predict(smooth_rate,
      #                           newdata =  data.frame(start = x[t]), se = TRUE)$fit
      x <- x + stepsize*predict(smooth_rate,
                                newdata =  data.frame(start = x), se = TRUE)$fit
    }
    
    return(x)
  }) %>% unlist
  
  otu_df <- otu_df %>% 
    mutate(predicted_mode = ifelse(end_prediction < change_point, "low", "high")) %>% 
    mutate(predicted_change = ifelse(start_mode != predicted_mode, "yes", "no")) %>% 
    mutate(otu = taxon)
  
  
  return(otu_df)
}


learn_drift <- function(delta_df, taxon = "Dialister", x_seq = NULL, grid_by = .01) {
  
  
  
  change_point <- change_points[taxon]
  
  # otu_abundance <- data.frame(abundance = t(abundances(atlas1006_log))[, taxon])
  
  otu_df <- delta_df %>% 
    filter(otu == taxon) %>% 
    arrange(start)
  
  # Fit polynomial to rates
  smooth_rate <- loess(unit_rate ~ start, otu_df)
  
  
  # Get predicted rate of change for each observation which isn't the subject's last observation
  if(is.null(x_seq)) {
    x_seq <- seq(from = range(otu_df$start)[1], to = range(otu_df$start)[2], by = grid_by)
  }
  
  
  rate_prediction <- predict(smooth_rate, newdata =  data.frame(start = x_seq), se = TRUE)
  
  rate_prediction_df <- data.frame(x = x_seq, 
                                   drift = rate_prediction$fit,
                                   se_drift = rate_prediction$se.fit)
  
  
  rate_prediction_df <- rate_prediction_df %>% 
    mutate(drift_lower = drift - se_drift, 
           drift_upper = drift + se_drift, 
           taxon = taxon)
  
  
  return(rate_prediction_df)
}

drift_residuals <- function(delta_df, taxon = "Dialister") {
  
  change_points <- c("Anaerofustis" = 0.4, 
                     "Aneurinibacillus" = 0.3, 
                     "Aquabacterium" = 0.5, 
                     "Bacteroides intestinalis et rel." = .5, 
                     "Bacteroides fragilis et rel." = 3.1,
                     "Burkholderia" = 0.5, 
                     "Dialister" = 3.25, 
                     "Leminorella" = 0.5, 
                     "Prevotella melaninogenica et rel." = 4.5, 
                     "Prevotella oralis et rel." = 4, 
                     "Serratia" = 0.5, 
                     "Uncultured Bacteroidetes" = 0.5, 
                     "Uncultured Clostridiales II" = 3.4, 
                     "Uncultured Selenomonadaceae" = 0.5, 
                     "Wissella et rel." = 0.45)
  
  change_point <- change_points[taxon]
  
  # otu_abundance <- data.frame(abundance = t(abundances(atlas1006_log))[, taxon])
  
  otu_df <- delta_df %>% 
    filter(otu == taxon) %>% 
    arrange(start)
  
  # Fit polynomial to rates
  smooth_rate <- loess(unit_rate ~ start, otu_df)
  
  
  res <- data.frame(residuals = smooth_rate$residuals,
                    start = otu_df$start, 
                    taxon = taxon)
  
  res <- res %>%
    mutate(abs_residuals = abs(residuals))
  
  
  return(res)
  
}

opper_diffusion <- function(delta_df, taxon = "Dialister", x_seq = NULL) {
  
  # Batz, Ruttor, Opper: Approximate Bayes learning of stochastic differential equations
  # p.5: diffusion can be approximated from (dX)^2/dt
  
  otu_df <- delta_df %>% 
    filter(otu == taxon) %>% 
    arrange(start) %>% 
    mutate(y_tilda = delta_abundance^2/delta_t)
  
  
  # Smooth
  smooth_diffusion <- loess(y_tilda ~ start, otu_df)
  
  
  # If no point is specified use a grid from first to last observations
  if(is.null(x_seq)) {
    x_seq <- seq(from = range(otu_df$start)[1], to = range(otu_df$start)[2], by = 0.01)
  }
  
  
  diffusion_prediction <- predict(smooth_diffusion, newdata =  data.frame(start = x_seq) %>% 
                                    set_colnames("start"), se = TRUE)
  
  diffusion_prediction_df <- data.frame(x = x_seq, 
                                        diffusion = diffusion_prediction$fit,
                                        se_diffusion = diffusion_prediction$se.fit)
  
  
  diffusion_prediction_df <- diffusion_prediction_df %>% 
    mutate(diffusion_lower = diffusion - se_diffusion, 
           diffusion_upper = diffusion + se_diffusion, 
           taxon = taxon)
  
  return(diffusion_prediction_df)
  
}
 

to_range <- function(x, range) {
  
  if(x < range[1]) {
    x <- range[1]
  } else if((x > range[2])) {
    x <- range[2] 
  } else {
    x <- x
  }
  
  return(x)
  
}

diffusion_mime <- function(data,
                           time = 1:90/30,
                           initial_value = NULL,
                           taxon = "Dialister",
                           detail = 10) {
  
  mydata <- data
  
  # res <- drift_residuals(mydata, taxon)
  drift <- learn_drift(mydata, taxon)
  
  taxon_range <- range(drift$x)
  
  
  # Initialize observation vector
  y <- rep(NA, length(time))
  
  # Random initial value. If out of range --> other end point
  if(is.null(initial_value)) {
    y[1] <- to_range(abs(rnorm(1, mean(taxon_range), 1)), range = taxon_range)
  } else y[1] <- to_range(initial_value, range = taxon_range)
  
  
  # Next values
  for(i in 2:length(time)) {
    
    # Divide each time step in to a finer grid to get more stable Euler approximation
    finer_x_grid <- seq(from = time[i-1], to = time[i], length.out = detail)
    
    temp_y <- y[i-1]
    
    # Finer grid values
    for(j in 2:length(finer_x_grid)) {
      
      # Time step
      dt <- finer_x_grid[j] - finer_x_grid[j-1]
      
      # Learn drift at this point
      point_drift <- learn_drift(mydata, taxon, temp_y)
      point_diffusion <- opper_diffusion(mydata, taxon, temp_y)
      
      
      # New value
      temp_y <- temp_y +
        point_drift$drift*dt +
        point_diffusion$diffusion*rnorm(1, 0, sqrt(dt))
      # point_drift$se_drift*(rt(1, 5)*sqrt(dt*3/5)) # Student version
      
      temp_y <- to_range(temp_y, taxon_range)
      
      
    }
    
    y[i] <- temp_y
    
  }
  
  return(data.frame(y = y, time = time))
}


diffusion_deterministic <- function(data,
                           time = 1:90/30,
                           initial_value = NULL,
                           taxon = "Dialister",
                           detail = 10) {
  
  mydata <- data
  
  # res <- drift_residuals(mydata, taxon)
  drift <- learn_drift(mydata, taxon)
  
  taxon_range <- range(drift$x)
  
  
  # Initialize observation vector
  y <- rep(NA, length(time))
  
  # Random initial value. If out of range --> other end point
  if(is.null(initial_value)) {
    y[1] <- to_range(abs(rnorm(1, mean(taxon_range), 1)), range = taxon_range)
  } else y[1] <- to_range(initial_value, range = taxon_range)
  
  
  # Next values
  for(i in 2:length(time)) {
    
    # Divide each time step in to a finer grid to get more stable Euler approximation
    finer_x_grid <- seq(from = time[i-1], to = time[i], length.out = detail)
    
    temp_y <- y[i-1]
    
    # Finer grid values
    for(j in 2:length(finer_x_grid)) {
      
      # Time step
      dt <- finer_x_grid[j] - finer_x_grid[j-1]
      
      # Learn drift at this point
      point_drift <- learn_drift(mydata, taxon, temp_y)
      # point_diffusion <- opper_diffusion(mydata, taxon, temp_y)
      point_diffusion <- 0 # Deterministic
      
      
      # New value
      temp_y <- temp_y +
        point_drift$drift*dt +
        point_diffusion*rnorm(1, 0, sqrt(dt))
      
      
      temp_y <- to_range(temp_y, taxon_range)
      
      
    }
    
    y[i] <- temp_y
    
  }
  
  return(data.frame(y = y, time = time))
}


range_seq <- function(x, increment = NULL, length_out = NULL) {
  
  if(is.null(increment) & is.null(length_out)) {
    stop("Choose increment or length_out (not both)")
  }
  
  # range
  r <- range(x)
  
  if(!is.null(increment)) {
    r_seq <- seq(from = r[1], to = r[2], by = increment)
  } else {
    r_seq <- seq(from = r[1], to = r[2], length.out = length_out)
  }
  
  
  return(r_seq)
  
}

## HitChip Atlas data ********************* ####

data("atlas1006")

# Log 10 transform
atlas1006_log <- atlas1006 %>% transform("compositional") %>% transform("log10")


# Subjects with multiple samples
multiple_id <- meta(atlas1006_log) %>% 
  filter(time != 0) %>% 
  pull(subject) %>% 
  unique()


longitudinal_atlas1006_log <- atlas1006_log %>% 
  subset_samples(subject %in% multiple_id)


longitudinal_atlas1006_log_abundances <- t(abundances(longitudinal_atlas1006_log))

## Bimodal taxa *************************** ####

# Bimodality coefficients
# atlas_1006_bimodality <- bimodality(longitudinal_atlas1006_log)

# Bimodal taxa
# bimodal_taxa <- names(which(atlas_1006_bimodality > 0.9))

bimodal_taxa <- c("Bacteroides fragilis et rel.", 
                  "Dialister", 
                  "Prevotella melaninogenica et rel.", 
                  "Prevotella oralis et rel.", 
                  "Uncultured Clostridiales II")



change_points <- c("Anaerofustis" = 0.4, 
                   "Aneurinibacillus" = 0.3, 
                   "Aquabacterium" = 0.5, 
                   "Bacteroides intestinalis et rel." = .5, 
                   "Bacteroides fragilis et rel." = 3.1,
                   "Burkholderia" = 0.5, 
                   "Dialister" = 3.25, 
                   "Leminorella" = 0.5, 
                   "Prevotella melaninogenica et rel." = 4.5, 
                   "Prevotella oralis et rel." = 4, 
                   "Serratia" = 0.5, 
                   "Uncultured Bacteroidetes" = 0.5, 
                   "Uncultured Clostridiales II" = 3.4, 
                   "Uncultured Selenomonadaceae" = 0.5, 
                   "Wissella et rel." = 0.45)

## Rates of change ************************ ####

# Get table with abundance changes 
delta_df <- lapply(meta(longitudinal_atlas1006_log) %>% pull(subject) %>% unique, function(i) {
  
  
  sample_info <- meta(longitudinal_atlas1006_log) %>% 
    filter(subject == i) %>% 
    select(time, sample)
  
  
  sample_abundances <- cbind(sample_info, longitudinal_atlas1006_log_abundances[sample_info$sample,
                                                             bimodal_taxa])
  
  # Average duplicated time index observations 
  duplicated_time <- sample_abundances$time[which(duplicated(sample_abundances$time))]
  if(length(duplicated_time) != 0) {
    for(d in duplicated_time) {
      mean_abundances <- sample_abundances[sample_abundances$time == d, bimodal_taxa] %>% 
        colMeans()
      
      sample_abundances[sample_abundances$time == d, bimodal_taxa][1, ] <- mean_abundances
      sample_abundances[sample_abundances$time == d, bimodal_taxa][2, ] <- NA
      
      sample_abundances <- sample_abundances %>% drop_na()
    }
  }
  
  
  deltas <- lapply(bimodal_taxa, function(otu) {
    otu_info <- sample_abundances[, c("time", otu)]
    
    delta <- lapply(1:(nrow(otu_info) - 1), function(r) {
      
       cbind(otu_info[r+1, ] - otu_info[r, ], start = otu_info[r, 2])
      
    }) %>% do.call(rbind, .)
     
    delta <- delta %>% 
      cbind(., otu) %>% 
      set_colnames(c("delta_t", "delta_abundance", "start", "otu"))
    
    return(delta)
    
  }) %>%
    do.call(rbind, .) %>% 
    cbind(subject = i)
  
  
  return(deltas)
  
}) %>% do.call(rbind, .)

# Add average change velocity
delta_df <- delta_df %>% 
  mutate(unit_rate = delta_abundance/delta_t)


otu_abundance <- data.frame(abundance = t(abundances(atlas1006_log))[, bimodal_taxa])

ggplot() + 
  geom_density(data = otu_abundance %>% 
                 melt() %>% 
                 set_colnames(c("sample", "otu", "abundance")), aes(x = abundance), fill = "red") + 
  geom_point(data = delta_df %>% 
               filter(otu %in% bimodal_taxa), aes(x = start, y = unit_rate), size = .75) + 
  geom_smooth(data = delta_df %>% 
                filter(otu %in% bimodal_taxa), aes(x = start, y = unit_rate)) +
  geom_hline(yintercept = 0, linetype = "dashed") + 
  facet_wrap(~otu, scales = "free")
  


## Predict future state with current 

# Each state which precedes an observations

# otu_ex <- "Dialister"
# change_point <- 3.25
# 
# otu_abundance <- data.frame(abundance = t(abundances(atlas1006_log))[, otu_ex])
# 
# otu_df <- delta_df %>% 
#   filter(otu == otu_ex) %>% 
#   mutate(end = start + delta_abundance) %>% 
#   mutate(start_mode = ifelse(start < change_point, "low", "high")) %>% 
#   mutate(end_mode = ifelse(start + delta_abundance < change_point, "low", "high")) %>% 
#   mutate(change = ifelse(start_mode != end_mode, "yes", "no"),
#          change_to_high = ifelse(start_mode == "low" & end_mode == "high", "yes", "no"), 
#          change_to_low = ifelse(start_mode == "high" & end_mode == "low", "yes", "no"))
# 
# 
# 
# ggplot() + 
#   geom_density(data = otu_abundance, aes(x = abundance), fill = "red") + 
#   geom_point(data = otu_df, aes(x = start, y = unit_rate)) + 
#   geom_smooth(data = otu_df, aes(x = start, y = unit_rate)) +
#   geom_hline(yintercept = 0, linetype = "dashed") + 
#   geom_vline(xintercept = 3.25, linetype = "dashed") +  
#   # geom_density(data = otu_df %>% filter(change_to_high == "yes"), aes(x = start), color = "green") +
#   # geom_density(data = otu_df %>% filter(change_to_low == "yes"), aes(x = start), color = "green") +
#   labs(x = "start", y = "Unit rate") +
#   facet_wrap(~otu) 
# 
# 
# 
# # Fit polynomial to rates
# smooth_rate <- loess(otu_df$unit_rate ~ otu_df$start)
# 
# 
# # Get predicted end point for each observation which isn't the subject's last observation
# end_prediction <- predict(smooth_rate, data.frame(start = otu_df$start), se = TRUE)
# 
# 
# otu_df$end_prediction <- otu_df$start + end_prediction$fit*otu_df$delta_t
# 
# otu_df <- otu_df %>% 
#   mutate(predicted_mode = ifelse(end_prediction < change_point, "low", "high")) %>% 
#   mutate(predicted_change = ifelse(start_mode != predicted_mode, "yes", "no"))
# 
# 
#  otu_df %>% 
#  ggplot() +
#   # geom_point(aes(x = end, y = end_prediction)) +
#   # geom_abline(intercept = 0, slope = 1)
#   geom_point(aes(x = start, y = end_prediction), color = "red") +
#   geom_point(aes(x = start, y = end))
#   
# 




## Inferred rates of change per taxa ****** ####







# Get results
premonitinon_res <- lapply(bimodal_taxa, function(taxon) {
  print(taxon)
  
  res <- bimodality_premonition(delta_df, taxon)
  
  return(res)
  
}) %>% do.call(rbind, .)


# Predicted mode changes
premonitinon_res <- premonitinon_res %>% 
  mutate(predicted_change_to_high = ifelse(start_mode == "low" & predicted_mode == "high", "yes", "no"), 
         predicted_change_to_low = ifelse(start_mode == "high" & predicted_mode == "low", "yes", "no"))
  


# Change in general
general_mode_accuracy <- premonitinon_res %>% 
  group_by(otu) %>% 
  summarize(mean(predicted_change == change, na.rm = T))

# Change from low to high
low_to_high_accuracy <- premonitinon_res %>% 
  group_by(otu) %>% 
  summarize(mean(predicted_change_to_high == change_to_high, na.rm = T))


# Change from high to_low
high_to_low_accuracy <- premonitinon_res %>% 
  group_by(otu) %>% 
  summarize(mean(predicted_change_to_low == change_to_low, na.rm = T))



genera_accuracy_df <- data.frame(general_mode_accuracy, low_to_high_accuracy[, 2], high_to_low_accuracy[, 2]) %>% 
  set_colnames(c("otu", "general", "high_to_low", "low_to_high")) %>% 
  arrange(general) %>% 
  mutate(otu = as.factor(otu))


# Mode prediction accuracy
genera_accuracy_df %>% 
  melt(id.vars = "otu") %>% 
  ggplot(aes(x = factor(otu, levels = rev(genera_accuracy_df$otu)), y = value, fill = variable)) +
  geom_col(position = "dodge") +
  theme_bw(15) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Genus", y = "Prediction accuracy")





premonitinon_res %>% 
  # filter(otu == bimodal_taxa[2]) %>% 
  ggplot() +
  geom_point(aes(x = start, y = end)) + 
  geom_point(aes(x = start, y = end_prediction), color = "red") + 
  facet_wrap(~otu, scales = "free")
  




  ggplot() +
  geom_density(data =  data.frame(t(abundances(atlas1006_log))[, bimodal_taxa]) %>% 
                 set_colnames(colnames(t(abundances(atlas1006_log))[, bimodal_taxa])) %>%
                 melt %>%
                 set_colnames(c("otu", "abundance")), aes(x = abundance), fill = "red") +
  geom_smooth(data = premonitinon_res %>%
               mutate(prediction_error = end - end_prediction), aes(x = start, y = abs(prediction_error))) +
    # geom_point(data = premonitinon_res, aes(x = start, y = 1)) +
    
  facet_wrap(~otu, scales = "free")
  
  
  
  
  
  premonitinon_res %>%
    mutate(prediction_error = end - end_prediction)%>% 
   group_by(otu) %>% 
    summarize(mean(prediction_error))
  
  

  # bimodal_taxa[-c(7, 9, 10, 16)]
  
ggplot() + 
  geom_density(data = t(abundances(atlas1006_log))[, bimodal_taxa] %>% 
                 melt() %>% 
                 as.data.frame() %>% 
                 set_colnames(c("sample", "otu", "value")), aes( x= value)) + 
  geom_histogram(data = premonitinon_res %>% 
                 filter(change == "yes"),
                 aes(x = start, y = ..count../max(..count..)),
                 color = "green", fill = "white", alpha = .1) + 
  facet_wrap(~otu, scales = "free")
  











## Learn drift **************************** ####



# Get drifts, residuals and diffusion for bimodal taxa
drift_res <- lapply(bimodal_taxa, function(taxon) {
  print(taxon)
  
  res <- learn_drift(delta_df, taxon)
  
  return(res)
  
}) %>% do.call(rbind, .)

drift_residual_res <- lapply(bimodal_taxa, function(taxon) {
  print(taxon)
  
  res <- drift_residuals(delta_df, taxon)
  
  return(res)
  
}) %>% do.call(rbind, .)

opper_res <- lapply(bimodal_taxa, function(taxon) {
  print(taxon)
  
  res <- opper_diffusion(delta_df, taxon)
  
  return(res)
  
}) %>% do.call(rbind, .)

# Potential functions
drift_res <- drift_res %>% 
  group_by(taxon) %>% 
  arrange(x) %>% 
  mutate(potential = cumsum(drift),
         potential_lower = cumsum(drift_lower),
         potential_upper = cumsum(drift_upper)) %>% 
  mutate(potential = (potential - min(potential))/max(potential)) %>%
  ungroup()


# Plot drift, averaged
p_drift <- ggplot() +
geom_hline(yintercept = 0, linetype = "dashed") +
  geom_density(data = otu_abundance %>%
                 melt %>%
                 set_colnames(c("sample", "taxon", "value")), 
               aes(x = value), fill = "grey", alpha = .5) +
geom_line(data = drift_res, aes(x = x, y = drift)) +
geom_ribbon(data = drift_res, aes(x = x, ymin = drift_lower, ymax = drift_upper), fill = "red", alpha = .25) +
  facet_wrap(~taxon, scales = "free", ncol = 1) +
labs(title = "Drift")

# 
# # Residuals of the drift fit
# p_volatility <- ggplot() +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   geom_density(data = otu_abundance %>%
#                  melt %>%
#                  set_colnames(c("sample", "taxon", "value")), 
#                aes(x = value), fill = "grey", alpha = .5) +
#   geom_point(data = drift_residual_res, aes(x = start, y = residuals), size = .25) +
#   geom_smooth(data = drift_residual_res, aes(x = start, y = (residuals)), size = .25, fill = "blue") +
#   facet_wrap(~taxon, scales = "free", ncol = 1) +
#   coord_cartesian(ylim = -1:1) +
#   labs(title = "Volatility")
# 
# 

# Opper diffusion
p_opper <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_density(data = otu_abundance %>%
                 melt %>%
                 set_colnames(c("sample", "taxon", "value")), 
               aes(x = value), fill = "grey", alpha = .5) +
  # geom_point(data = opper_res, aes(x = x, y = diffusion), size = .25) +
  geom_line(data = opper_res, aes(x = x, y = diffusion), size = .25) +
  geom_ribbon(data = opper_res, aes(x = x, ymin = diffusion_lower, ymax = diffusion_upper),
              fill = "green", alpha = .25) +
  facet_wrap(~taxon, scales = "free", ncol = 1) +
  coord_cartesian(ylim = -1:1) +
  labs(title = "Diffusion")


plot_grid(p_drift, p_opper,  ncol = 2)
  
  

## Diffusion mime ************************* ####




taxon <- bimodal_taxa[2]

# mime_seires <- diffusion_mime(delta_df,
#                               time = 1:90/30,
#                               initial_value = 3,
#                               taxon = taxon,
#                               detail = 10)


time_span <- (1:(360*10))/30 

mime_set_of_sets <- lapply(5, function(otu) {
  
  otu_range <- otu_abundance[, otu] %>% range
  d_range <- otu_range[2] - otu_range[1]
  
  mime_set <- lapply(otu_range[2], function(i) {
    
    mime_seires <- diffusion_mime(delta_df,
                                  time = time_span,
                                  initial_value = i,
                                  taxon = taxon,
                                  detail = 5)
    
    
    return(mime_seires$y)
  }) %>% 
    do.call(cbind, .) %>% 
    cbind(., time = time_span, taxon = otu) %>% 
    data.frame()
  
  
  return(mime_set)
  
}) %>% do.call(rbind, .)

mime_set_of_sets <- mime_set_of_sets %>% 
  data.frame()

 

mime_set_of_sets[, 1] <- mime_set_of_sets[, 1] %>% as.character %>% as.numeric
mime_set_of_sets[, 2] <- mime_set_of_sets[, 2] %>% as.character %>% as.numeric
mime_set_of_sets[, 3] <- mime_set_of_sets[, 3] %>% as.character %>% as.numeric
mime_set_of_sets[, 4] <- mime_set_of_sets[, 4] %>% as.character %>% as.numeric


  mime_set_of_sets %>% 
    filter(taxon == 5, time <= 100) %>% 
    filter((time%% .1) == 0) %>% 
    # drop_na() %>% 
  melt(id.vars = c("taxon", "time")) %>%
  # filter(taxon == "Dialister") %>% 
  ggplot(aes(x = time, y = value)) +
  geom_line(size = 1) +
  # geom_hline(data = as.data.frame(change_points) %>%
  #              rownames_to_column(var = "taxon") %>%
  #              filter(taxon %in% bimodal_taxa),
  #            aes(yintercept = change_points), linetype = "dashed") +
  # facet_wrap(~taxon, ncol = 1, scales = "free") +
    guides(color = FALSE) + 
    labs(title ="Simulation", x = "Time", y = "")
  
q
  
## Early warning signals ****************** ####

mime_set_of_sets %>% 
    filter(taxon == bimodal_taxa[3]) %>% 
    select(-c(V2, V3)) %>% 
    ggplot(aes(x = time, y = V1)) +
    geom_line(color = "red") +
    geom_hline(yintercept = change_points[bimodal_taxa[3]], linetype = "dashed") + 
    geom_vline(xintercept = 9.1, linetype = "dashed", color = "red") +
    labs(subtitle = "Example series with change point")
  
  
  mime_set_of_sets %>% 
    filter(taxon == bimodal_taxa[3], 
           time < 9.1) %>% 
    pull(V1) %>%
    generic_ews(detrending = "gaussian")
  
  
  
  
delta_df %>%
  ggplot(aes(x = start, y = 1/unit_rate)) +
  geom_smooth() +
  facet_wrap(~otu, scales = "free")


## Empirical acf(1) *********************** ####

taxon <- bimodal_taxa[2]
acf_at <- 1

acf1 <- lapply(range_seq(otu_abundance[, taxon], increment = .1), function(x) {
    
  # Deterministic value
  x_deterministic <- diffusion_deterministic(data = delta_df,
                                             time = 0:1,
                                             initial_value = x,
                                             taxon = taxon,
                                             detail = 10)
  
  x_deterministic <- x_deterministic[x_deterministic$time == 1, "y"]
  
  
  # Simulate n_realizations realizations from a given x
  n_realizations <- 50
  x_diffusion <- lapply(1:n_realizations, function(i) {
    print(i)
    sim_val <- diffusion_mime(data = delta_df,
                   time = 0:1,
                   initial_value = x,
                   taxon = taxon,
                   detail = 10)
    
    
    
    sim_val[sim_val$time == acf_at, "y"]
    
  }) %>%
    unlist %>% 
    data.frame %>% 
    cbind(y = ., index = 1:n_realizations, x = x, type = "stochastic") %>% 
    set_colnames(c("y", "index", "x", "type"))
    
  
  
  
  # Empirical acf(1)
  data.frame(x = x, ac = mean((x_diffusion$y - x_deterministic)^2))
  
  
}) %>% do.call(rbind, .)


dialister_dispersion <- opper_diffusion(delta_df, taxon, range_seq(otu_abundance[, taxon], increment = .1))




xx <- full_join(dialister_dispersion, acf1)

xx <- xx %>% mutate(ac =ac/diffusion)
  
  ggplot() + 
  geom_density(data = otu_abundance[, taxon] %>% data.frame(abundance = .), 
               aes(x = abundance)) +
  geom_point(data = xx, aes(x = x, y = ac), color = "black", size = 1.5) + 
    geom_point(data = xx, aes(x = x, y = ac), color = "red", size = .5) + 
  geom_smooth(data = xx,aes(x = x, y = ac), color = "red") +
    geom_smooth(data = dialister_dispersion, aes(x = x, y = diffusion)) +
    labs(title = taxon, subtitle = "blue = dispersion, red = simulated acf(1)")
  

  
  
  xx <- xx %>% mutate(inv_ac = 1/ac) %>%
    mutate(inv_ac = inv_ac/max(inv_ac, na.rm = T))
  
  ggplot() + 
    geom_density(data = otu_abundance[, taxon] %>% data.frame(abundance = .), 
                 aes(x = abundance)) +
    # geom_point(data = xx, aes(x = x, y = ac), color = "black", size = 1.5) + 
    # geom_point(data = xx, aes(x = x, y = ac), color = "red", size = .5) + 
    geom_smooth(data = xx,aes(x = x, y = inv_ac), color = "red") +
    geom_smooth(data = dialister_dispersion, aes(x = x, y = diffusion)) +
    labs(title = taxon, subtitle = "blue = dispersion, red = simulated acf(1)")
  
  
## Loop drift and dispersion ************** ####

all_taxa <- top_taxa(longitudinal_atlas1006_log, n = 20)
all_abundance <- t(abundances(atlas1006_log))[, all_taxa]

all_delta_df <- lapply(meta(longitudinal_atlas1006_log) %>% pull(subject) %>% unique, function(i) {
  
  
  sample_info <- meta(longitudinal_atlas1006_log) %>% 
    filter(subject == i) %>% 
    select(time, sample)
  
  
  sample_abundances <- cbind(sample_info, longitudinal_atlas1006_log_abundances[sample_info$sample,
                                                                                all_taxa])
  
  # Average duplicated time index observations 
  duplicated_time <- sample_abundances$time[which(duplicated(sample_abundances$time))]
  if(length(duplicated_time) != 0) {
    for(d in duplicated_time) {
      mean_abundances <- sample_abundances[sample_abundances$time == d, all_taxa] %>% 
        colMeans()
      
      sample_abundances[sample_abundances$time == d, all_taxa][1, ] <- mean_abundances
      sample_abundances[sample_abundances$time == d, all_taxa][2, ] <- NA
      
      sample_abundances <- sample_abundances %>% drop_na()
    }
  }
  
  
  deltas <- lapply(all_taxa, function(otu) {
    otu_info <- sample_abundances[, c("time", otu)]
    
    delta <- lapply(1:(nrow(otu_info) - 1), function(r) {
      
      cbind(otu_info[r+1, ] - otu_info[r, ], start = otu_info[r, 2])
      
    }) %>% do.call(rbind, .)
    
    delta <- delta %>% 
      cbind(., otu) %>% 
      set_colnames(c("delta_t", "delta_abundance", "start", "otu"))
    
    return(delta)
    
  }) %>%
    do.call(rbind, .) %>% 
    cbind(subject = i)
  
  
  return(deltas)
  
}) %>% do.call(rbind, .)

# Add average change velocity
all_delta_df <- all_delta_df %>% 
  mutate(unit_rate = delta_abundance/delta_t)

  

# Get drifts, residuals and diffusion for bimodal taxa
all_drift_res <- lapply(all_taxa, function(taxon) {
  print(taxon)
  
  res <- learn_drift(all_delta_df, taxon)
  
  return(res)
  
}) %>% do.call(rbind, .)

all_opper_res <- lapply(all_taxa, function(taxon) {
  print(taxon)
  
  res <- opper_diffusion(all_delta_df, taxon)
  
  return(res)
  
}) %>% do.call(rbind, .)




# Potential functions
all_drift_res <- all_drift_res %>% 
  group_by(taxon) %>% 
  arrange(x) %>% 
  mutate(potential = cumsum(drift),
         potential_lower = cumsum(drift_lower),
         potential_upper = cumsum(drift_upper)) %>% 
  mutate(potential = (potential - min(potential))/max(potential)) %>%
  ungroup()


# Plot drift, averaged
p_drift <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_density(data = all_abundance %>%
                 melt %>%
                 set_colnames(c("sample", "taxon", "value")), 
               aes(x = value), fill = "grey", alpha = .5) +
  geom_line(data = all_drift_res, aes(x = x, y = drift)) +
  geom_ribbon(data = all_drift_res, aes(x = x, ymin = drift_lower, ymax = drift_upper), fill = "red", alpha = .25) +
  facet_wrap(~taxon, scales = "free") +
  labs(title = "Drift")


# Opper diffusion
p_opper <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_density(data = all_abundance %>%
                 melt %>%
                 set_colnames(c("sample", "taxon", "value")), 
               aes(x = value), fill = "grey", alpha = .5) +
  # geom_point(data = opper_res, aes(x = x, y = diffusion), size = .25) +
  geom_line(data = all_opper_res, aes(x = x, y = diffusion*x), size = .25) +
  geom_ribbon(data = all_opper_res, aes(x = x, ymin = diffusion_lower, ymax = diffusion_upper),
              fill = "green", alpha = .25) +
  facet_wrap(~taxon, scales = "free") +
  coord_cartesian(ylim = -1:1) +
  labs(title = "Diffusion")


plot_grid(p_drift, p_opper,  ncol = 2)
## Nonparametric stationary density ******* ####

otu <- all_taxa[1]

drift_diff_df <- full_join(all_drift_res %>% 
  filter(taxon == otu) %>% 
  select(x, drift), 
all_opper_res %>% 
  filter(taxon == otu) %>% 
  select(x, diffusion), 
  by = "x")

nonparametric_stationry_density <- function(drift, diffusion, x_seq) {
  
  x0 <- x_seq[1]
  
  # Log scale measure
  integrand <- drift/(diffusion^2)
  log_s <- -2*cumsum(integrand)
  
  # Speed measure
  log_m <- - log(diffusion^2) - log_s
  
  
  stationary_density <- exp(log_m)
  
  # Normalize
  stationary_density <- stationary_density/sum(stationary_density)
  
  return(data.frame(density = stationary_density, x = x_seq))
}




dyn_density <- nonparametric_stationry_density(drift_diff_df$drift, drift_diff_df$diffusion, drift_diff_df$x)


dyn_density

ggplot() + 
  geom_density(data = data.frame(abundance = all_abundance[, otu]),
               aes(x = abundance)) +
  geom_line(data = dyn_density, aes(x = x, y = density*100), color = "red")
