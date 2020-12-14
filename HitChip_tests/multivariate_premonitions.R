## Functions  ***************************** ####

# Function to get all unique pairs of vector elements
get_pairs <- function(x) {
  
  save_class <- class(x)
  
  unique_x <- unique(x) %>% as.character
  
  x_pairs <- list()
  
  row_ind <- 1
  
  for(i in 1:(length(unique_x) - 1)) {
    for(j in (i+1):length(unique_x)) {
      
      x_pairs[[row_ind]] <- c(unique_x[i], unique_x[j])
      
      row_ind <- row_ind + 1
      
    }
  }
  
  return(x_pairs)
} 

# Function to get pairwise starting points and abundances
get_pairwise_df <- function(delta_df, otu1, otu2) {
  
  x_start <- delta_df %>% 
    filter(otu == otu1) %>% 
    pull(start)
  
  x_change <- delta_df %>% 
    filter(otu == otu1) %>% 
    pull(unit_rate)
  
  
  y_start <- delta_df %>% 
    filter(otu == otu2) %>% 
    pull(start)
  
  y_change <- delta_df %>% 
    filter(otu == otu2) %>% 
    pull(unit_rate)
  
  
  pairwise_df <- data.frame(x_start, x_change, y_start, y_change, id = 1:length(x_start))
  
  return(pairwise_df)
  
}


# Set up grid, essentially a wrapper for expand.grid
make_grid <- function(x_range, y_range,
                      resolution = 10, padding = 0) {
  
  
  x_grid <- seq(from = x_range[1] - padding, to = x_range[2] + padding, length.out = resolution) 
  y_grid <- seq(from = y_range[1] - padding, to = y_range[2] + padding, length.out = resolution)
  
  
  df <- data.frame(x = rep(x_grid, resolution), 
                   y = rep(y_grid, each = resolution))
  
  return(df)
}

# Get grid averges for a single variable z
grid_average <- function(x, y, z, resolution = 10, padding = 0) {
  
  mydf <- data.frame(x, y, z)
  
  mydf <- mydf %>% mutate(index = 1:nrow(mydf))
  
  # Get xy ranges
  x_range <- x %>% range
  y_range <- y %>% range
  xy_grid <- make_grid(x_range, y_range, resolution, padding)
  
  # Add column for cell averages
  xy_grid$z <- NA
  
  
  for(i in 1:nrow(xy_grid)) {
    
    # cell observations
    cell_df <- mydf %>% 
      filter(x < xy_grid[i, 1] & y < xy_grid[i, 2])
    
    # cell average
    xy_grid[i, "z"] <- mean(cell_df$z)
    
    # Remove used observations from data
    mydf <- mydf %>% filter(!(index %in% cell_df$index))
  }
  
  # Replace NaN --> NA
  xy_grid[is.nan(xy_grid$z), "z"] <- NA
  
  
  return(xy_grid)
}


multiple_2d_loess <- function(x, y, z, newdata = NULL, resolution = 10, span = .1) {
  
  z <- data.frame(z)
  
  
  # If no new data is specified, use even grid in the x,y-range
  if(is.null(newdata)) {
    x_range <- range(x)
    y_range <- range(y)
    
    newdata <- expand.grid(x = seq(from = x_range[1], to = x_range[2], length.out = resolution), 
                           y = seq(from = y_range[1], to = y_range[2], length.out = resolution))
  }
  
  # Loop over z components (= vector components)
  zs <- lapply(1:ncol(z), function(k) {
    # Fit loess
    fit <- loess(z[, k] ~ x + y, span)
    
    # Predict on new data.
    prediction <- predict(fit, newdata, se = TRUE)
    
    # Get predicted values
    predicted_z <- prediction$fit %>%
      melt %>%
      pull(value)
    
    # Get predicted standard errors
    predicted_se <- prediction$se %>%
      melt %>%
      pull(value)
    
    
    return(data.frame(predicted_z, predicted_se))
  }) %>%
    do.call(cbind, .) %>% 
    as.data.frame() %>% 
    set_colnames(paste0(rep(colnames(z), each = ncol(z)), c("", "_se")))
  

  # Results
  res <- cbind(newdata, zs)
  
  return(res)
}


to_range2d <- function(vec, x_range, y_range) {
  
  vec[1] <- to_range(vec[1], x_range)
  vec[2] <- to_range(vec[2], y_range)
  
  return(vec)
  
}

## HitChip Atlas data ********************* ####

data("atlas1006")

# Log 10 transform
atlas1006_log <- atlas1006 %>% transform("log10")


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
    
    # Add time step ID
    
    delta <- delta %>% 
      mutate(time_id = 1:nrow(delta))
    
    return(delta)
    
  }) %>%
    do.call(rbind, .) %>% 
    cbind(subject = i)
  
  
  return(deltas)
  
}) %>% do.call(rbind, .)

# Add average change velocity
delta_df <- delta_df %>% 
  mutate(unit_rate = delta_abundance/delta_t)

pairwise_df <- get_pairwise_df(delta_df, bimodal_taxa[1], bimodal_taxa[2])

# Plot observations
p <- pairwise_df %>% 
  ggplot() +
  geom_point(aes(x = x_start, y = y_start), color = "red", size = 1) +
  geom_segment(aes(x = x_start, y = y_start,
                   xend = x_start + x_change, yend = y_start + y_change, group = id),
               arrow = arrow(length = unit(0.2, "cm"), ends="last"), size = .5)





## What is this *************************** ####
gl <- 40

x_fit <- interp.loess(pairwise_df$x_start, pairwise_df$y_start, pairwise_df$x_change, gridlen = c(gl, gl))
y_fit <- interp.loess(pairwise_df$x_start, pairwise_df$y_start, pairwise_df$y_change, gridlen = c(gl, gl))

# data.frame(x_start = x_fit$x, y_start = x_fit$y, 
#            xend = x_fit$z, yend = y_fit$z) %>% head

df <- data.frame(x_change = (data.frame(x_fit$z) %>% melt)[, 2], 
y_change = (data.frame(y_fit$z) %>% melt)[, 2], 
x = rep(x_fit$x, gl), 
y = rep(x_fit$y, each = gl))


q <- df %>% 
  ggplot() +
  geom_point(aes(x = x, y = y), color = "red", size = .25) +
  geom_segment(aes(x = x, y =y,
                   xend = x + x_change/50, yend = y + y_change/50),
               arrow = arrow(length = unit(0.1, "cm"), ends="last"), size = .5) +
  coord_cartesian(xlim = 2.5:3.5)




plot_grid(p, q)





## Grid averages ************************** ####

## Empirical 

p_empirical_changes <- ggplot() +
  geom_density_2d(data = otu_abundance %>%
                    as.data.frame() %>%
                    set_colnames(paste0("V", 1:5)), aes(x = V1, y = V2), color = "grey") +
  geom_segment(data = pairwise_df, aes(x = x_start, y = y_start,
                                  xend =  x_start + x_change/5, yend = y_start + y_change/5),
               arrow = arrow(length = unit(0.15, "cm"), ends="last"), size = .5) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red") + 
  labs(subtitle = "Empirical changes")+ 
  coord_cartesian(ylim = 2:5)


##### Average rate in a grid





mean_x <- grid_average(pairwise_df$x_start, pairwise_df$y_start, pairwise_df$x_change, resolution = 15, padding = 0)

mean_y <- grid_average(pairwise_df$x_start, pairwise_df$y_start, pairwise_df$y_change, resolution = 15, padding = 0)


mean_xy <- full_join(mean_y, mean_x, by = c("x", "y")) %>% 
  set_colnames(c("x", "y", "x_change", "y_change"))


# Plot 
p_emricial_grid_avegare <-   ggplot() + 
  geom_density_2d(data = otu_abundance %>%
                    as.data.frame() %>%
                    set_colnames(paste0("V", 1:5)), aes(x = V1, y = V2), color = "grey") +
  geom_segment(data =mean_xy, aes(x = x, y = y,
                   xend =  x + x_change/5, yend = y + y_change/5),
               arrow = arrow(length = unit(0.15, "cm"), ends="last"), size = .5) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red") + 
  labs(subtitle = "Average change per time unit, 15x15 grid")+ 
  coord_cartesian(ylim = 2:5)


taxon1_drift <- learn_drift(delta_df, taxon = bimodal_taxa[1], grid_by = .1)
taxon2_drift <- learn_drift(delta_df, taxon = bimodal_taxa[2], grid_by = .1)


# taxon1_drift[, c("x", "drift")] %>% set_colnames(c("x", "x_drift"))
# taxon2_drift[, c("x", "drift")] %>% set_colnames(c("y", "y_drift"))

xy_grid <- expand.grid(x = taxon1_drift[, "x"], y = taxon2_drift[, "x"])
xy_grid$x_drift <- NA
xy_grid$y_drift <- NA

# Get drifts
for(i in 1:nrow(xy_grid)) {
  print(i)
  xy_grid[i, "x_drift"] <- taxon1_drift %>% filter(x == xy_grid[i, "x"]) %>% pull(drift)
  xy_grid[i, "y_drift"] <- taxon2_drift %>% filter(x == xy_grid[i, "y"]) %>% pull(drift)
  
}


p_independent_drifts_learned <- ggplot() +
  geom_density_2d(data = otu_abundance %>%
                    as.data.frame() %>%
                    set_colnames(paste0("V", 1:5)), aes(x = V1, y = V2), color = "grey") +
  geom_segment(data = xy_grid, aes(x = x, y = y,
                                  xend =  x + x_drift/3, yend = y + y_drift/3),
               arrow = arrow(length = unit(0.1, "cm"), ends="last"), size = .5) + 
  labs(subtitle = "Individual learned drifts") + 
  coord_cartesian(ylim = 2:5)



plot_grid(p_empirical_changes, p_emricial_grid_avegare, p_independent_drifts_learned, p_2d_drifts_learned, nrow = 1)

## Interpolate 1D valued function ********* #### 




x_fit <- interp.loess(mean_xy$x, mean_xy$y, mean_xy$x_change, gridlen = c(gl, gl))
y_fit <- interp.loess(mean_xy$x, mean_xy$y, mean_xy$y_change, gridlen = c(gl, gl))

df <- data.frame(x_change = (data.frame(x_fit$z) %>% melt)[, 2], 
                 y_change = (data.frame(y_fit$z) %>% melt)[, 2], 
                 x = rep(x_fit$x, gl), 
                 y = rep(x_fit$y, each = gl))



  ggplot() +
  geom_density_2d(data = otu_abundance %>% 
                    as.data.frame() %>% 
                    set_colnames(paste0("V", 1:5)), aes(x = V2, y = V3)) + 
  geom_segment(data = df, aes(x = x, y =y,
                   xend = x + x_change/50, yend = y + y_change/50),
               arrow = arrow(length = unit(0.1, "cm"), ends="last"), size = .5)

  
  
  ggplot() +
    geom_density_2d(data = otu_abundance %>% 
                      as.data.frame() %>% 
                      set_colnames(paste0("V", 1:5)), aes(x = V2, y = V3)) + 
    geom_segment(data = mean_xy, aes(x = x, y =y,
                                xend = x + x_change, yend = y + y_change),
                 arrow = arrow(length = unit(0.1, "cm"), ends="last"), size = .5)

  
  

## 2D loess drift ************************* ####
  
# Fit loess and predict on new data




mfit <- multiple_2d_loess(x = pairwise_df$x_start, y = pairwise_df$y_start, z = pairwise_df[, c(2,4)],
                          resolution = 30, span = .5, 
                          newdata = make_grid(range(otu_abundance[, bimodal_taxa[1]]), 
                                              range(otu_abundance[, bimodal_taxa[2]]), resolution = 20))


# Combine standard errors. This is probably not an ok way
# mfit <- mfit %>%
#   mutate(se = x_change_se + y_change_se)
  

# PLOT
p_2d_drifts_learned <- ggplot() + 
geom_density_2d(data = otu_abundance %>%
                  as.data.frame() %>%
                  set_colnames(paste0("V", 1:5)), aes(x = V1, y = V2), color = "grey") +
  # geom_point(data = otu_abundance %>%
  #                   as.data.frame() %>%
  #                    set_colnames(paste0("V", 1:5)), aes(x = V1, y = V2), size = .2) +
  # geom_point(data = mfit, aes(x = x, y= y), color = "red", size = .2) +
  geom_segment(data = mfit, aes(x = x, y =y,
                            xend = x + x_change/2.5, yend = y + y_change/2.5),
              arrow = arrow(length = unit(0.1, "cm"), ends="last"), size = .5) +
  labs(subtitle = "2D learned drifts") + 
  coord_cartesian(ylim = 2:5)




## 2D loess dispersion ******************** ####


diffusion_2d_loess <- function(x, y, z, newdata = NULL, resolution = 10, span = .1) {
  
  z <- data.frame(z)
  
  # Square differencies 
  z <- apply(z, 2, FUN = function(i) (i^2))
  
  # If no new data is specified, use even grid in the x,y-range
  if(is.null(newdata)) {
    x_range <- range(x)
    y_range <- range(y)
    
    newdata <- expand.grid(x = seq(from = x_range[1], to = x_range[2], length.out = resolution), 
                           y = seq(from = y_range[1], to = y_range[2], length.out = resolution))
  }
  
  # Loop over z components (= vector components)
  zs <- lapply(1:ncol(z), function(k) {
    # Fit loess
    fit <- loess(z[, k] ~ x + y, span)
    
    # Predict on new data.
    prediction <- predict(fit, newdata, se = TRUE)
    
    # Get predicted values
    predicted_z <- prediction$fit %>%
      melt %>%
      pull(value)
    
    ## Replace negative values with 0
    predicted_z[predicted_z < 0] <- 0
    
    # Get predicted standard errors
    predicted_se <- prediction$se %>%
      melt %>%
      pull(value)
    
    
    return(data.frame(predicted_z, predicted_se))
  }) %>%
    do.call(cbind, .) %>% 
    as.data.frame() %>% 
    set_colnames(c("x_diffusion", "x_diffusion_se", "y_diffusion", "y_diffusion_se"))
  
  
  # Results
  res <- cbind(newdata, zs)
  
  return(res)
}

ggplot() + 
  geom_density_2d(data = otu_abundance %>%
                    as.data.frame() %>%
                    set_colnames(paste0("V", 1:5)), aes(x = V1, y = V2)) +
  # geom_point(data = otu_abundance %>%
  #                   as.data.frame() %>%
  #                    set_colnames(paste0("V", 1:5)), aes(x = V1, y = V2), size = .2) +
  # geom_point(data = mfit, aes(x = x, y= y), color = "red", size = .2) +
  geom_point(data = dispersion, aes(x = x, y =y))



## Diffusion mime 2 *********************** ####

diffusion_mime_2d <- function(data,
                           time = 1:90/30,
                           initial_value = NULL,
                           taxon1 = "Bacteroides fragilis et rel.",
                           taxon2 = "Dialister",
                           detail = 10, 
                           resolution = 10,
                           span = .5) {
  
  # mydata <- delta_df
  mydata <- data
  
  pairwise_df <- get_pairwise_df(mydata, taxon1, taxon2)
  
  #####
  # For some reason this doesn't respect the x,y ranges. Could be an extrapolation issue
  ##### 
  drift_fit <- multiple_2d_loess(x = pairwise_df$x_start, y = pairwise_df$y_start,
                            z = pairwise_df[, c("x_change", "y_change")],
                            resolution, span, 
                            newdata = make_grid(range(otu_abundance[, taxon1]), 
                                                range(otu_abundance[, taxon2]), resolution))
  
  # Remove NA 
  drift_fit <- drift_fit %>% drop_na
  
  taxon_range <- data.frame(x = range(drift_fit$x), y = range(drift_fit$y))
  
  
  # Initialize observation matrix
  obs <- matrix(NA,  length(time), 2)
  
  # Random initial value. If out of range --> other end point
  if(is.null(initial_value)) {
    obs[1, ] <- to_range2d(c(abs(rnorm(1, mean(taxon_range$x), 1)), 
                           abs(rnorm(1, mean(taxon_range$y), 1))),
                           taxon_range$x, taxon_range$y)
  } else {
    obs[1, ] <- to_range2d(initial_value, taxon_range$x, taxon_range$y)
  } 
    
  
  
  # Next values *************************
  for(i in 2:length(time)) {
    print(i)
    # Divide each time step in to a finer grid to get more stable Euler approximation
    finer_x_grid <- seq(from = time[i-1], to = time[i], length.out = detail)
    
    temp_y <- obs[i-1, ]
    
    # Finer grid values
    for(j in 2:length(finer_x_grid)) {
    
      
      # Time step
      dt <- finer_x_grid[j] - finer_x_grid[j-1]
      
      # Learn drift at this point
      point_drift <- multiple_2d_loess(x = pairwise_df$x_start, y = pairwise_df$y_start,
                                       z = pairwise_df[, c(2,4)],
                                       resolution = resolution, span = span, 
                                       newdata = data.frame(x = temp_y[1], y = temp_y[2]) %>%
                                         set_colnames(c("x", "y")))
      
      point_diffusion <- diffusion_2d_loess(x = pairwise_df$x_start, y = pairwise_df$y_start,
                                          z = pairwise_df[, c(2,4)],
                                          resolution = resolution, span = span, 
                                          newdata = data.frame(x = temp_y[1], y = temp_y[2]) %>%
                                            set_colnames(c("x", "y")))
      
      
      # New value
      temp_y <- temp_y +
        point_drift[, c("x_change", "y_change")]*dt +
        point_diffusion[, c("x_diffusion", "y_diffusion")]*rnorm(1, 0, sqrt(dt))
      
      temp_y <- to_range2d(temp_y, taxon_range$x, taxon_range$y)
      
      
    }
    
    obs[i, ] <- temp_y %>% as.numeric
    
  }
  
  res <- data.frame(x = obs[, 1], y = obs[, 2], time = time)
  
  return(res)
}



xx <- diffusion_mime_2d(delta_df, time = 1:90/30, initial_value = c(3.5, 2.5), resolution = 20, span = .1)

ggplot() +
  geom_path(data =  xx, aes(x = x, y = y, color = time))

ggplot() +
  geom_path(data =  xx, aes(x = time, y= y))
ggplot() +
  geom_path(data =  xx, aes(x = time, y= x))


ggplot() + 
  geom_density_2d(data = otu_abundance %>%
                    as.data.frame() %>%
                    set_colnames(paste0("V", 1:5)), aes(x = V1, y = V2), color = "grey") +
  geom_path(data =  xx, aes(x = x, y = y, color = time*30)) +
  scale_color_gradient(low = "blue", high = "red") +
  theme_classic(15) +
  labs(x = "Bacteroides fragilis et rel.", y = "Dialister", title = "Simulation") +
  guides(color = guide_legend(title = "Time"))
  






diffusion_fit <- diffusion_2d_loess(x = pairwise_df$x_start, y = pairwise_df$y_start,
                               z = pairwise_df[, c("x_change", "y_change")],
                               resolution, span = 1.5, 
                               newdata = make_grid(range(otu_abundance[, taxon1]), 
                                                   range(otu_abundance[, taxon2]), resolution = 50))

