# HitChip functions

## Packages
library(tidyverse)
library(microbiome)
library(magrittr)
library(reshape2)

## General *********************************** ####



# Return x to range
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

# Return x to 2D range
to_range2d <- function(vec, x_range, y_range) {
  
  vec[1] <- to_range(vec[1], x_range)
  vec[2] <- to_range(vec[2], y_range)
  
  return(vec)
  
}

# Get a sequence along the range of x
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
get_pairwise_df <- function(delta_df, x_otu, y_otu) {
  
  x <- delta_df %>% 
    filter(otu == x_otu) %>% 
    select(x, unit_dx)
  
  y <- delta_df %>% 
    filter(otu == y_otu) %>% 
    select(x, unit_dx)
  
  
  
  pairwise_df <- data.frame(x[, "x"], 
                            y[, "x"],
                            x[, "unit_dx"],
                            y[, "unit_dx"]) %>% 
    set_colnames(c("x", "y", "unit_dx", "unit_dy"))
  
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


# Corner of a data frame
corner <- function(x, n = 5) {
  
  x[1:n, 1:n]
  
}


## Data getters ****************************** ####

# These functions read and prepare the data for gaussian process interpolators

# There are dublicated time indices in hitchip data. This functions selects the sample with higher read count and returns the edited phyloseq
get_edit_hitchip_pseq <- function(pseq) {
  
  # data("atlas1006")
  # pseq <- atlas1006
  
  # Edit taxa_names
  atlas_taxa <- taxa_names(pseq) %>% 
    gsub(" ", "_", .)
  taxa_names(pseq) <- atlas_taxa
  
  
  # Add read counts
  sample_data(pseq) <- cbind(meta(pseq), read_count = rowSums(t(abundances(pseq))))
  
  # Subjects with multiple samples
  multiple_id <- meta(pseq) %>% 
    filter(time != 0) %>% 
    pull(subject) %>% 
    unique()
  
  # Mark all longitudinal samples as dublicates, some of these are removed
  sample_data(pseq) <- cbind(meta(pseq),
                                  dub = meta(pseq) %>% 
    mutate(dub = ifelse(subject %in% multiple_id, "yes", "no")) %>% 
    pull(dub))
  
  meta <- meta(pseq)
  
  for(i in unique(meta$subject)) {
    
    # Time indices
    subject_time <- meta %>% 
      filter(subject == i) %>% 
      pull(time)
    
    # Dublicatd times  
    dub_time <- subject_time[which(duplicated(subject_time))]
    
    # Remove smaller read counts from meta
    for(t in unique(dub_time)) {
      
      max_read_count <- meta %>% 
        filter(subject == i, time == t) %>% 
        pull(read_count) %>% 
        max()
      
      meta[meta$subject == i &
             meta$time == t &
             (meta$read_count != max_read_count), ] <- NA
      
      meta <- meta %>% drop_na
      
    }
    
    
  }
  
  
  sample_data(pseq) <- meta
  
  return(pseq)
}

get_hitchip_data <- function(pseq, transformation = "log10", pcoa = FALSE) {
  

  # Edit taxa_names
  atlas_taxa <- taxa_names(pseq)
  
  
  # Subjects with multiple samples
  multiple_id <- meta(pseq) %>% 
    filter(time != 0) %>% 
    pull(subject) %>% 
    unique()
  
  
  
  # Log 10 transform
  pseq_log <- pseq %>%
    microbiome::transform("compositional") %>%
    microbiome::transform(transformation)
  
  # Subset pseq to longitudinal samples
  longitudinal_pseq_log <- pseq_log %>% 
    subset_samples(subject %in% multiple_id)
  
  # Longitudinal log abundances
  longitudinal_pseq_log_abundances <- t(abundances(longitudinal_pseq_log)) %>% 
    set_colnames(atlas_taxa)
  
  
  
  
  # Get table with abundance changes 
  print("Taxon changes")
  delta_df <- lapply(meta(longitudinal_pseq_log) %>%
                       pull(subject) %>%
                       unique, function(i) {
                         
                         
                         sample_info <- meta(longitudinal_pseq_log) %>% 
                           filter(subject == i) %>% 
                           select(time, sample)
                         
                         
                         sample_abundances <- cbind(sample_info, longitudinal_pseq_log_abundances[sample_info$sample, ])
                         
                         # # Average duplicated time index observations 
                         # duplicated_time <- sample_abundances$time[which(duplicated(sample_abundances$time))]
                         # if(length(duplicated_time) != 0) {
                         #   for(d in duplicated_time) {
                         #     mean_abundances <- sample_abundances[sample_abundances$time == d, atlas_taxa] %>% 
                         #       colMeans()
                         #     
                         #     sample_abundances[sample_abundances$time == d, atlas_taxa][1, ] <- mean_abundances
                         #     sample_abundances[sample_abundances$time == d, atlas_taxa][2, ] <- NA
                         #     
                         #     sample_abundances <- sample_abundances %>% drop_na()
                         #   }
                         # }
                         
                         
                         deltas <- lapply(atlas_taxa, function(otu) {
                           otu_info <- sample_abundances[, c("time", otu)]
                           
                           delta <- lapply(1:(nrow(otu_info) - 1), function(r) {
                             
                             cbind(otu_info[r+1, ] - otu_info[r, ], start = otu_info[r, 2])
                             
                           }) %>% do.call(rbind, .)
                           
                           delta <- delta %>% 
                             cbind(., otu) %>% 
                             set_colnames(c("dt", "dx", "x", "otu"))
                           
                           return(delta)
                           
                         }) %>%
                           do.call(rbind, .) %>% 
                           cbind(subject = i)
                         
                         
                         return(deltas)
                         
                       }) %>% 
    do.call(rbind, .)
  
  
  # Add unit velocity
  delta_df <- delta_df %>% 
    mutate(unit_dx = dx/dt)
  
  
  full_abundance_table <- data.frame(t(abundances(pseq_log)))
  
  # initialize results list
  return_list <- list()
  
  # Add PCoA
  if(pcoa) {
    print("PCoA")
    ord <- ordinate(pseq, "MDS", "bray")
    
    # Take first 2 PCoA aexs
    pcoa_df <- ord$vectors[, 1:2]
    
    # Longitudinal PCoA samples; add identifiers as well
    longitudinal_pcoa_df <- cbind(pcoa_df[sample_names(longitudinal_pseq_log), ], 
                                  meta(longitudinal_pseq_log)[, c("time", "subject")])
    
    pcoa_delta_df <- lapply(unique(longitudinal_pcoa_df$subject), function(i) {
      
      sample_abundances <- longitudinal_pcoa_df %>% 
        filter(subject == i)
      
      
      deltas <- lapply(c("Axis.1", "Axis.2"), function(axis) {
        
        otu_info <- sample_abundances[, c("time", axis)]
        
        
        delta <- lapply(1:(nrow(otu_info) - 1), function(r) {
          
          cbind(otu_info[r+1, ] - otu_info[r, ], start = otu_info[r, 2])
          
        }) %>% do.call(rbind, .)
        
        
        delta <- delta %>% 
          cbind(., axis) %>% 
          set_colnames(c("dt", "dx", "x", "axis"))
        
        return(delta)
        
      }) %>%
        do.call(rbind, .) %>% 
        cbind(subject = i) %>% 
        mutate(unit_dx = dx /dt)
      
      
      return(deltas)
      
    }) %>% 
      do.call(rbind, .)
    
    
    return_list[["pcoa"]] <- list(delta_df = pcoa_delta_df, 
                                  full_pcoa = pcoa_df)
  }
  
  
  return_list[["delta_df"]] <- delta_df
  return_list[["abundances"]] <- full_abundance_table
  return_list[["longitudinal_abundances"]] <- longitudinal_pseq_log_abundances
  
  return(return_list)
  
}








get_david_data <- function(transformation = "log10") {
  
  # Read data
  load("/Users/villelaitinen/Desktop/PhD/early_warning_signals/data-David2014/David_phyloseq.Rdata")
  
  david_names <- c("A", "B")

  # Transform pseqs
  pseq_transformed <- lapply(list(seqA, seqB), function(pseq) {
    
    # Transform
    trans_pseq <- pseq %>%
      aggregate_taxa("Genus") %>%
      microbiome::transform("compositional") %>%
      microbiome::transform(transformation)
    
    return(trans_pseq)
  }) %>%
    set_names(david_names)
  
  david_taxa <- taxa_names(pseq_transformed[[1]])
  
  # Get abundances
  otus <- lapply(pseq_transformed, function(pseq) {
    
    # Get abundances and add time
    abn <- t(abundances(pseq)) %>% 
      as.data.frame %>% 
      mutate(time = meta(pseq)$time)
    
    return(abn)
  }) %>% 
    set_names(david_names)
 
  
  # Get abundance differencies
  deltas <- lapply(otus, function(otu_tbl) {
    
    # Otu differences
    dotu <- otu_tbl[2:nrow(otu_tbl), ] - otu_tbl[1:(nrow(otu_tbl) - 1), ] 
    
    # Change time column name
    dotu <- dotu %>% 
      mutate(dt = time) %>% 
      select(-time)
    
    return(dotu)
  }) %>%
    set_names(david_names)
  
  # Edit to standard form
  delta_df <- lapply(david_names, function(id) {
    
    edit_delta <- deltas[[id]] %>% 
      melt(id.vars = "dt") %>% 
      transmute(dt = dt,
                otu = variable,
                dx = value, 
                x = otus[[id]][1:(nrow(otus[[id]]) - 1), ] %>%
                  melt(id.vars = "time") %>% 
                  pull(value), 
                unit_dx = dx/dt) 
    
    return(edit_delta)
    
  }) %>%
    set_names(david_names)
  
  
  return(list(delta_df = delta_df, 
              adundances = do.call(rbind, otus) %>% 
                select(-time)))
  
}
get_ravel_data <- function() {
  
  load(file = "/Users/villelaitinen/Desktop/PhD/early_warning_signals/data-Ravel2012/data.clr.rda")
  
  N_subjects <- dim(data.cube.clr)[1]
  
  otu_df <- lapply(1:N_subjects, function(i) {
    df <- data.cube.clr[i, , ]
    
    # Remove white space from column names
    colnames(df) <- gsub(" ", "", colnames(df))
    
    # Add time (rownames) to column
    
    df <- cbind(df, time =  rownames(df) %>% as.numeric)
    
    return(df)
  })
  
  deltas <- lapply(otu_df, function(df) {
    
    # Abundance differencies, change time column to dt
    delta_abundances <- as.data.frame(df[2:nrow(df), ] - df[1:(nrow(df) - 1), ]) %>% 
      mutate(dt = time) %>% 
      select(-time)
    
    # Melt
    delta_abundances <- delta_abundances %>% 
      melt(id.vars = "dt") %>% 
      transmute(dt = dt, 
                otu = variable, 
                dx = value)
    
    return(delta_abundances)
    
  })
  
  
  delta_df <- lapply(1:N_subjects, function(i) {
    
    start_x_and_time <- otu_df[[i]][1:(nrow(otu_df[[i]]) - 1), ] %>% 
      as.data.frame %>%
      melt(id.vars = "time") %>% 
      select(time, value)
    
    
    df <- deltas[[i]] %>% 
      mutate(x = start_x_and_time$value,
             start_time = start_x_and_time$time,
             unit_dx = dx / dt, 
             subject = i)
    
    
    df <- df %>% 
      group_by(otu) %>% 
      mutate(rounded_dx = round(dx, 10)) %>% 
      mutate(duplicate = duplicated(rounded_dx)) %>% 
        select(-rounded_dx) %>% 
      ungroup
    
    
    
    return(df)
    
  }) %>% 
    do.call(rbind, .)
  
  # Times at which a 
  duplicate_start_times <- lapply(1:N_subjects, function(i) {
    delta_df %>% 
      filter(subject == i, duplicate == TRUE) %>% 
      pull(start_time) %>% 
      unique()
  })
  
  # Remove duplicates: duplicate_start_times gives the indeces to remove
  # delta_df
  no_duplicates_delta_df <- lapply(1:N_subjects, function(i) {
  
    # Remove imputed values
    unique_deltas <- delta_df %>% 
      filter(subject == i) %>% 
      filter(!(start_time %in% duplicate_start_times[[i]])) %>% 
      select(-duplicate)
    
    return(unique_deltas)
    
  }) %>% 
    do.call(rbind, .)
  # abundance table
  no_duplicates_otu_df <- lapply(1:N_subjects, function(i) {
    
    # Remove imputed values
    unique_otus <- otu_df[[i]] %>% 
      data.frame() %>% 
      filter(!(time %in% duplicate_start_times[[i]])) %>% 
      select(-time)
      
    
    return(unique_otus)
    
  })
  
  
  
  return(list(delta_df = no_duplicates_delta_df, 
              abundances = do.call(rbind, no_duplicates_otu_df)
              )
  )
  
  
  
  
  # 
  # pca <- prcomp(ravel_df_combined)
  # 
  # pca_clr_time <- cbind(pca$x[, 1:2],
  #                       time =  as.numeric(ravel_df_combined %>% rownames),
  #                       subject = rep(1:32, each = 95)) %>% 
  #   as.data.frame()
  # 
  # rownames(pca_clr_time) <- NULL

  
}


# ravel_data <- get_ravel_data()
# pairwise_df <- get_pairwise_df(ravel_data$delta_df, "Dialister", "Aerococcus")

## Loess fitter ****************************** ####

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

# Get estimate for drift term
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

# Get estimates for diffusion term
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



# Simulate diffusion process
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

# Simulate deterministic counterpart
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





# 2D loess interpolator
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



## GP imputation ***************************** ####

# Implement vector values gaussian process interpolator as per
# [1] Alvarez et al: Kernels for vector-valued functions: a Review 


# Standard square exponential kernel
se_covariance <- function(x1, x2, sigma, rho, add_linear = FALSE, mean_c, linear_factor = .5) {
  
  N1 <- length(x1)
  N2 <- length(x2)
  
  res <- matrix(NA, N1, N2)
  
  for(i in 1:N1) {
    for(j in 1:N2) {
      delta_x_sq <- (x1[i] - x2[j])^2
      
      
      
      res[i, j] <- sigma*exp(-0.5*sum( delta_x_sq / rho^2 ))
    }
  }
  
  
  if(add_linear)  {
    
    for(i in 1:N1) {
      for(j in 1:N2) {
        
        res[i, j] <- res[i, j] + linear_factor*(x1[i] - mean_c)*(x2[j] - mean_c)
          
      }
    }
    
  }
  
  
  
  return(res)
}

# Use decomposable kernel as in [2] Yildiz et al.: Learning stochastic differential equations with gaussian processes without gradient matching 
sub_K <- function(x1, x2, sigma, rho, D, mean_c, add_linear = FALSE) {
  
  # Stop if input vector not in same space
  if(length(x1) != length(x2)) {
    stop("Lengths of input vectors do not match!")
  }
  
  
  # Squared exponential
  delta_x_sq <- (x1 - x2)^2
  res <- sigma*exp(-0.5*sum( delta_x_sq / rho^2 ) )
  
  # OUP
  # delta_x <- (x1 - x2)
  # res <- sigma*exp(sum( delta_x / rho ) )

  # Add linear kernel to avoid patological extrapolation
  if(add_linear)  {
    lin_res <- lapply(1:length(x1), function(i) {
      (x1[i] - mean_c)*(x2[i] - mean_c)
    }) %>% unlist
    
    res <-  res + 0.5*lin_res
  }
  

  return(res*diag(D))
  
}

# Helper: builds block matrix from list of blocks
block_list_to_matrix <- function(block_list) {
  
  # First index level denotes the row number
  res <- lapply(block_list, function(i) {
    do.call(cbind, i)
  }) %>% do.call(rbind, .)
  
  return(res)
} 

# Build block matrix. X1 is matrix, where rows are the observation locations and X2 are prediction locations
K <- function(X1, X2, sigma, rho, D, mean_c, add_linear = FALSE) {
  
  if(class(X1) == "numeric") {
    X1 <- t(as.matrix(X1))
  }
  if(class(X2) == "numeric") {
    X2 <- t(as.matrix(X2))
  }
  
  N1 <- nrow(X1) # number of X1 points
  N2 <- nrow(X2) # number of X2 points
  
  ## Make list of block matrices
  # First list level is for rows
  block_list <- lapply(1:N1, function(i) {
    lapply(1:N2, function(j) {
      sub_K(X1[i, ], X2[j, ], sigma, rho, D, mean_c, add_linear)
    })
  })
  
  
  block_matrix <- block_list_to_matrix(block_list)
  
  return(block_matrix)
  
}

# Make matrix for observation noise, with standard deviation = epsilon
noise_matrix <- function(epsilon, N) {
  
  D <- length(epsilon)
  M_sigma <- epsilon*diag(D)
  
  kronecker_sigma <- M_sigma %x% diag(N)
  
  return(kronecker_sigma)
  
}

## Predictive mean and covariance.
# X = observation locations, Y = observation values, x_pred = prediction location, 
# epsilon = noise standard deviation
f_pred <- function(X, Y, x_pred, sigma, rho, epsilon) {
  
  N <- nrow(X)
  
  
  # K1 <- K(X, X, sigma, rho); print("K1") # This is time consuming
  K2 <- K(x_pred, X, sigma, rho) 
  K3 <- K(X, x_pred, sigma, rho) 
  K4 <- K(x_pred, x_pred, sigma, rho) 
  noise <- noise_matrix(epsilon, N)
  
  # Build formula 8 in [1]
  
  inv_K_plus_noise <- solve(K1 + noise)
  
  pred_mean <- (K2 %*% inv_K_plus_noise) %*% as.vector(as.matrix(Y))
  pred_mean <- pred_mean %>% as.vector
  pred_covariance <- K4 - K2 %*% inv_K_plus_noise %*% K3
  
  return(list(mean = pred_mean, covariance = pred_covariance))
}




GP_imputation <- function(X, Y, X_pred, 
                          rho, sigma, epsilon, 
                          return_X_covariance = FALSE, 
                          K1 = NULL, 
                          mean_c = 0,
                          add_linear = FALSE, 
                          n_samples = 20, 
                          WAIC_samples = 1000) {
  
  # To data frame
  if(class(X) == "numeric") X <- data.frame(X)
  if(class(Y) == "numeric") Y <- data.frame(Y)
  if(class(X_pred) == "numeric") X_pred <- data.frame(X_pred)
  
  # Input & output dimension
  D_input <- ncol(X)
  D_output <- ncol(Y)
  
  
  # Number of training points
  N <- nrow(X)
  
  ## Check GP parameters ****************************
  if(length(mean_c) == 1) mean_c <- rep(mean_c, D_input)
  if(length(rho) == 1) rho <- rep(rho, D_input)
  # if(length(sigma) == 1) sigma <- rep(sigma, D_input)
  if(length(epsilon) == 1) epsilon <- rep(epsilon, D_output)
  
  if(length(rho) != D_input) stop("Invalid length-scale")
  if(length(sigma) != 1) stop("Invalid GP variance")
  if(length(epsilon) != D_output) stop("Invalid observations noise")
  
  
  
  
  
  ## Predictions ************************************
  # Build formula (8) in [1]
  K1 <- K(X, X, sigma, rho, D_output, mean_c, add_linear)
  print("K1") # This is time consuming but can be optimized further
  
  noise <- noise_matrix(epsilon^2, N)
  inv_K_plus_noise <- solve(K1 + noise)
  
  # Predictions: Normal distributions
  pred_list <- lapply(1:nrow(X_pred), function(i) {
    
    pred_location <- X_pred[i, ]
    
    K2 <- K(pred_location, X, sigma, rho, D_output, mean_c, add_linear) 
    K3 <- K(X, pred_location, sigma, rho, D_output, mean_c, add_linear) 
    K4 <- K(pred_location, pred_location, sigma, rho, D_output, mean_c, add_linear) 
    
    pred_mean <- (K2 %*% inv_K_plus_noise) %*% as.vector(as.matrix(Y))
    pred_mean <- pred_mean %>% as.vector
    pred_covariance <- K4 - K2 %*% inv_K_plus_noise %*% K3
    
    return(list(mean = pred_mean, covariance = pred_covariance))
    
  })
  # Prediction means
  pred_mean <- lapply(pred_list, function(i) {
    
    i$mean
    
  }) %>% do.call(rbind, .)
  # Prediction covariances
  pred_covariance <- lapply(pred_list, function(i) {
    
    i$covariance
    
  }) 
  
  
  ## WAIC *******************************************
  train_mean <- se_covariance(X[, 1], X[, 1], sigma, rho,
                              add_linear, mean_c,
                              linear_factor = .5) %*%
    inv_K_plus_noise %*%
    as.matrix(Y)
  
  sigma_star <- se_covariance(X[, 1], X[, 1], sigma, rho,
                              add_linear, mean_c,
                              linear_factor = .5) -
    t(se_covariance(X[, 1], X[, 1], sigma, rho, 
                    add_linear, mean_c,
                    linear_factor = .5)) %*%
    inv_K_plus_noise %*% 
    se_covariance(X[, 1], X[, 1], sigma, rho, 
                  add_linear, mean_c,
                  linear_factor = .5)
  
  
  # Samples
  function_samples <- rmvnorm(n = WAIC_samples, 
                              mean = train_mean[, 1], 
                              sigma = sigma_star)
  # Add noise
  # function_samples <- function_samples + rnorm(WAIC_samples*ncol(function_samples), 0, epsilon)
  
  # point-wise likelihoods
  ll <- apply(function_samples, 1, FUN = function(i) dnorm(Y[, 1], i, epsilon)) %>% t
  
  # Add nugget to combat underflow
  ll[ll == 0] <- 1e-12
  
  # WAIC
  lppd <- sum(log(colMeans(ll)))
  p_WAIC <- sum(apply(ll, 2, FUN = function(i) var(log(i))))
  
  WAIC <- -2*(lppd - p_WAIC)
  
  ## Function samples ******************************
  
  # Prediction covariance
  sigma_star <- se_covariance(X_pred[, 1], X_pred[, 1], sigma, rho,
                              add_linear, mean_c,
                              linear_factor = .5) -
    t(se_covariance(X[, 1], X_pred[, 1], sigma, rho, 
                    add_linear, mean_c,
                    linear_factor = .5)) %*%
    inv_K_plus_noise %*% 
    se_covariance(X[, 1], X_pred[, 1], sigma, rho, 
                  add_linear, mean_c,
                  linear_factor = .5)

  # Samples
  function_samples <- rmvnorm(n = n_samples, 
                              mean = pred_mean[, 1], 
                              sigma = sigma_star)
  
  

  
  
  ## Build results **********************************
  res <- list(prediction = list(mean = pred_mean, 
                                covariance = pred_covariance, 
                                X = X_pred), 
              training = list(X = X, Y = Y), 
              parameters = list(rho = rho, sigma = sigma, epsilon = epsilon), 
              samples = list(X = X_pred, samples = function_samples), 
              WAIC = WAIC)
  
  
  # Computation of K1 takes a lot of time
  # so might be useful to return it as well for consecutive predictions
  if(return_X_covariance) res$X_covariance <- K1
  
  return(res)
  
  
  
}


GP_2D_diffusion_plotter <- function(GP_drift, GP_diffusion, full_data, arrow_length = .15, arrow_scaling, variables) {
  
  full_data <- full_data %>% 
    data.frame %>% set_colnames(c("x", "y"))
  
  ## Training data
  training_df <- cbind(GP_drift$training$X, GP_drift$training$Y) %>% 
    set_colnames(c("x", "y", "unit_dx", "unit_dy"))
  
  
  
  # Base
  p_base <- ggplot() +
    geom_point(data = full_data, aes(x = x, y = y), size = .5, color = "blue") +
    geom_point(data = training_df, aes(x = x, y = y), size = .75) +
    geom_segment(data = training_df, aes(x = x, y = y, xend = x + unit_dx, yend = y + unit_dy), 
                 arrow = arrow(length = unit(arrow_length, "cm"), ends="last"), size = .25)
    
    
  
  
  ## Drift 
  drift_df <- cbind(GP_drift$prediction$X, GP_drift$prediction$mean) %>% 
    data.frame %>% 
    set_colnames(c("x", "y", "unit_dx", "unit_dy")) %>% 
    mutate(unit_dx = unit_dx*arrow_scaling, 
           unit_dy = unit_dy*arrow_scaling)
  
  # Add standard deviation
  drift_df$standard_deviation <- lapply(GP_drift$prediction$covariance, function(i) {
    sqrt(i[1, 1])
  }) %>% unlist

  
  # plot
  p_drift <-ggplot() +
    geom_segment(data = drift_df, aes(x = x, y = y, xend = x + unit_dx, yend = y + unit_dy),
                 arrow = arrow(length = unit(arrow_length, "cm"), ends="last"), size = .25)
  
  
  
  # Drift uncertainty
  p_drift_variance <- ggplot() +
    geom_raster(data = drift_df, aes(x = x, y = y, fill = standard_deviation)) +
    geom_density_2d(data = full_data, aes(x = x, y = y), size = .5, color = "blue") +
    # geom_point(data = full_data, aes(x = x, y = y), size = .5, color = "blue") +
    scale_fill_continuous(low = "white", high = "black", name = "Drift variance") + 
    labs(title = "Uncertainty on drift")
    
  
  # ## Base and drift
  # p_base_drift <- ggplot() +
  #   geom_point(data = training_df, aes(x = x, y = y), size = .75, color = "grey") +
  #   geom_segment(data = training_df, aes(x = x, y = y, xend = x + unit_dx, yend = y + unit_dy), 
  #                arrow = arrow(length = unit(0.15, "cm"), ends="last"), size = .25, color = "grey") +
  #   geom_point(data = drift_df, aes(x = x, y = y), size = .75, color = "red") +
  #   geom_segment(data = drift_df, aes(x = x, y = y, xend = x + unit_dx, yend = y + unit_dy), 
  #                arrow = arrow(length = unit(arrow_length, "cm"), ends="last"), size = .4)
  
  
  ## Diffusion
  
  diffusion_df <- data.frame(X = GP_diffusion$prediction$X,
                             mean = GP_diffusion$prediction$mean,
                        covariance = do.call(rbind, GP_diffusion$prediction$covariance) %>% sqrt)%>% 
    set_colnames(c("x", "y", "mean", "standard_deviation"))
                 
  p_diffusion <- ggplot() + 
    geom_raster(data = diffusion_df, aes(x = x, y = y, fill = mean)) +
    geom_density_2d(data = full_data, aes(x = x, y = y), size = .5, color = "blue") +
    scale_fill_continuous(low = "white", high = "red", name = "Diffusion")+ 
    labs(title = "Diffusion magnitude")
  
  # Uncertainty of diffusion
  p_diffusion_variance <- ggplot() + 
    geom_raster(data = diffusion_df, aes(x = x, y = y, fill = standard_deviation)) +
    geom_density_2d(data = full_data, aes(x = x, y = y), size = .5, color = "blue") +
    scale_color_continuous(low = "blue", high = "red") +
    scale_fill_continuous(low = "white", high = "black", name = "Diffusion variance") + 
    labs(title = "Diffusion uncertainty")
  
  
  ## Drift and diffusion
  p_drift_diffusion <- ggplot() + 
    geom_raster(data = diffusion_df, aes(x = x, y = y, fill = mean)) +
    # geom_density_2d(data = full_data, aes(x = x, y = y), size = .5, color = "blue") +
    geom_point(data = full_data, aes(x = x, y = y), size = .5, color = "blue") +
    geom_segment(data = drift_df, aes(x = x, y = y, xend = x + unit_dx, yend = y + unit_dy), 
                 arrow = arrow(length = unit(arrow_length, "cm"), ends="last"), size = .25) +
    scale_fill_continuous(low = "white", high = "red", name = "Diffusion")+ 
    labs(title = "Drift and Diffusion")
  
  
  return(list(base = p_base, 
              drift = p_drift, 
              diffusion = p_diffusion, 
              diffusion_variance = p_diffusion_variance, 
              drift_diffusion = p_drift_diffusion, 
              drift_variance = p_drift_variance))
}


GP_1D_diffusion_plotter <- function(GP_drift, GP_diffusion, full_data, variable, scaling = 1) {
  
  full_data <- full_data %>% 
    data.frame %>% set_colnames(c("x"))
  
  ## Training data
  training_df <- cbind(GP_drift$training$X, GP_drift$training$Y) %>% 
    set_colnames(c("x", "unit_dx")) %>% 
    mutate(unit_dx = unit_dx*scaling)
  
  ## Base ************************************************************
  p_base <- ggplot() +
    geom_point(data = training_df, aes(x = x, y = unit_dx), size = .75) +
    geom_density(data = full_data, aes(x = x), size = .5) +
     labs(subtitle = variable)
    
  
  
  
  
  ## Drift ************************************************************

  drift_df <- cbind(GP_drift$prediction$X, GP_drift$prediction$mean) %>%
    data.frame %>%
    set_colnames(c("x", "mean")) %>% 
    mutate(mean = mean*scaling)
  
  # Drift standard deviation
  drift_df$standard_deviation <- lapply(GP_drift$prediction$covariance, function(i) {
    sqrt(i[1, 1]*scaling)
  }) %>% unlist

  p_drift <- ggplot() +
    geom_density(data = full_data, aes(x = x), size = .5) +
    geom_ribbon(data = drift_df, aes(x = x, ymin = mean - standard_deviation, ymax = mean + standard_deviation), fill = "blue", alpha = .5) +
    geom_line(data = drift_df, aes(x = x, y = mean), size = .75) +
    labs(title = "Drift", subtitle = variable)


  ## Diffusion ************************************************************
  
  diffusion_df <- data.frame(X = (GP_diffusion$prediction$X),
                             mean = GP_diffusion$prediction$mean,
                             covariance = sqrt(do.call(rbind, GP_diffusion$prediction$covariance)*scaling)) %>% 
    set_colnames(c("x", "mean", "standard_deviation")) %>% 
    mutate(mean = mean*scaling)
  
  p_diffusion <- ggplot() +
    geom_density(data = full_data, aes(x = x), size = .5) +
    geom_ribbon(data = diffusion_df, aes(x = x, ymin = mean - standard_deviation, ymax = mean + standard_deviation), fill = "red", alpha = .5) +
    geom_line(data = diffusion_df, aes(x = x, y = mean), size = .75) +
    labs(title = "Diffusion", subtitle = variable)
  


  ## Drift and diffusion ************************************************************
  diffusion_drift_df <- rbind(cbind(diffusion_df, type = "Diffusion"), cbind(drift_df, type = "Drift"))
  
  p_drift_diffusion <- ggplot() + 
    geom_density(data = full_data, aes(x = x), size = .5) +
    geom_ribbon(data = diffusion_drift_df,
                aes(x = x, ymin = mean - standard_deviation,
                    ymax = mean + standard_deviation, 
                    fill = type), alpha = .5) +
    geom_line(data = diffusion_drift_df,
              aes(x = x, y = mean, group = type)) +
    labs(title = "Drift and Diffusion", subtitle = variable)
  
  
  return(list(base = p_base, 
              drift = p_drift,
              diffusion = p_diffusion, 
              drift_diffusion = p_drift_diffusion))
  
  
}


GP_1D_diffusion_plotter_set <- function(res_list, full_data, variable, scaling = 1) {
  
  full_data <- full_data %>% 
    data.frame %>% set_colnames(c("x"))
  
  # Combine drift tables 
  drift_df <- lapply(1:length(res_list), function(i) {
    
    GP_drift <- res_list[[i]]$drift
    GP_diffusion <- res_list[[i]]$diffusion
    
      
      drift_df <- cbind(GP_drift$prediction$X, GP_drift$prediction$mean) %>%
      data.frame %>%
      set_colnames(c("x", "mean")) %>% 
      mutate(mean = mean*scaling)
    
    drift_df$standard_deviation <- lapply(GP_drift$prediction$covariance, function(i) {
      sqrt(i[1, 1]*scaling)
    }) %>% unlist
      
    return(cbind(drift_df, subject = i))
  }) %>%
    do.call(rbind, .)
  
  # Combine diffusion tables 
  diffusion_df <- lapply(1:length(res_list), function(i) {
    
    GP_drift <- res_list[[i]]$drift
    GP_diffusion <- res_list[[i]]$diffusion
    
    
    diffusion_df <- data.frame(X = (GP_diffusion$prediction$X),
                               mean = GP_diffusion$prediction$mean,
                               covariance = do.call(rbind, GP_diffusion$prediction$covariance)*scaling %>% sqrt()) %>% 
      set_colnames(c("x", "mean", "standard_deviation")) %>% 
      mutate(mean = mean*scaling)
    
    
    return(cbind(diffusion_df, subject = i))
  }) %>%
    do.call(rbind, .)
  
  
  
  # Make drift plot
  p_drift <- ggplot() + 
    geom_density(data = full_data, aes(x = x), size = .5) +
    geom_ribbon(data = drift_df,
                aes(x = x, ymin = mean - standard_deviation,
                    ymax = mean + standard_deviation, 
                    group = subject), fill = "blue",  alpha = .1) +
    geom_line(data = drift_df,
              aes(x = x, y = mean, group = subject)) +
    labs(title = "Drift and Diffusion", subtitle = variable)
  
  # Make diffusion plot
  p_diffusion <- ggplot() + 
    geom_density(data = full_data, aes(x = x), size = .5) +
    geom_ribbon(data = diffusion_df,
                aes(x = x, ymin = mean - standard_deviation,
                    ymax = mean + standard_deviation, 
                    group = subject), fill = "red",  alpha = .1) +
    geom_line(data = diffusion_df,
              aes(x = x, y = mean, group = subject)) +
    labs(title = "Drift and Diffusion", subtitle = variable)
  
  
  return(list(drift = p_drift, 
              diffusion = p_diffusion))
  
}

add_difference_sq_2D <- function(pairwise_df) {
  
  
  res <- pairwise_df %>% 
    mutate(x_end = x + unit_dx, 
           y_end = y + unit_dy) %>%
    mutate(difference_sq = (x-x_end)^2 + (y-y_end)^2)

    
  return(res)
}
add_difference_sq_1D <- function(delta_df) {
  
  mydf <- delta_df
  
  mydf <- mydf %>% 
    mutate(difference_sq = dx^2)
  
  return(mydf)
  
}





hitchip_stationary_density_plotter <- function(hitchip_data, hitchip_otu,
                                               drift_pars, diffusion_pars) {
  
  
  drift_rho <- 1
  drift_sigma <- .1
  drift_epsilon <- .7
  
  diffusion_rho <- 1.5
  diffusion_sigma <- 1
  diffusion_epsilon <- 0.3
  
  
  
  x_pred <- range_seq(range(hitchip_data$abundances[, hitchip_otu])[1:2], length_out = 200)
  
  mean_c <- mean(hitchip_data$abundances[, hitchip_otu])
  
  n_samples <- 50
  
  ## One dimensional
  data <- hitchip_data$delta_df %>% 
    filter(otu == hitchip_otu) %>% 
    filter(dx != 0) 
  
  one_dimensional_hitchip_drift <- GP_imputation(X = data %>% 
                                                   pull(x), 
                                                 Y = data %>% 
                                                   pull(unit_dx),
                                                 X_pred = x_pred,
                                                 rho = drift_rho,
                                                 sigma = drift_sigma,
                                                 epsilon = drift_epsilon, 
                                                 mean_c = mean_c,
                                                 add_linear = FALSE, 
                                                 n_samples = n_samples)
  
  one_dimensional_hitchip_diffusion <- GP_imputation(X = data %>%
                                                       pull(x),
                                                     Y = data %>%
                                                       pull(dx) %>%
                                                       abs,
                                                     X_pred = x_pred,
                                                     rho = diffusion_rho,
                                                     sigma = diffusion_sigma,
                                                     epsilon = diffusion_epsilon,
                                                     mean_c, add_linear = FALSE,
                                                     n_samples = n_samples)
  
  
  
  
  hitchip_1d_plots <- GP_1D_diffusion_plotter(GP_drift = one_dimensional_hitchip_drift,
                                              GP_diffusion = one_dimensional_hitchip_diffusion,
                                              full_data = hitchip_data$abundances[, hitchip_otu],
                                              variable = hitchip_otu)
  
  
  
  ## Stationary densities
  
  # Average
  drift <- one_dimensional_hitchip_drift$prediction$mean
  diffusion <- one_dimensional_hitchip_diffusion$prediction$mean
  mean_density <- SDE_non_stationary_density(x_pred, drift, diffusion)
  
  # Samples
  sample_densities <- lapply(1:n_samples, function(i) {
    
    X <- one_dimensional_hitchip_drift$samples$X[, 1]
    
    drift <- one_dimensional_hitchip_drift$samples$samples[i, ]
    diffusion <- one_dimensional_hitchip_diffusion$samples$samples[i, ]
    
    res <- SDE_non_stationary_density(X, drift, diffusion)
    
    res <- res %>% 
      mutate(sample = i)
    
    return(res)
  }) %>% do.call(rbind, .)
  
  
  
  p <- ggplot() + 
    stat_density(data = data.frame(x = hitchip_data$abundances[, hitchip_otu]), 
                 aes(x = x), geom = "line") +
    geom_line(data = sample_densities, 
              aes(x = x, y = y, group = sample), color = "orange", alpha = 0.25) +
    geom_line(data = mean_density, 
              aes(x = x, y = y), color = "brown1") +
    labs(subtitle = hitchip_otu)
  
  
  return(p)
  
}



# p_zoom_drift <- loop_plot$drift + labs(subtitle = "Ravel; Bifidobacterium drift; all subjects; zoomed") + coord_cartesian(xlim = -1:0.5, ylim = -1:3)
# 
# p_big_drift <- loop_plot$drift + labs(subtitle = "Ravel; Bifidobacterium drift; all subjects")
# 
# p_zoom_diffusion <- loop_plot$diffusion + labs(subtitle = "Ravel; Bifidobacterium diffusion; all subjects; zoomed") + coord_cartesian(xlim = -1:0.5, ylim = -1:3)
# 
# p_big_diffusion <- loop_plot$diffusion + labs(subtitle = "Ravel; Bifidobacterium diffusion; all subjects")
# 
# 
# plot_grid(p_zoom_drift, p_big_drift, p_zoom_diffusion, p_big_diffusion, ncol = 2)


# 
# sigma <- .1
# rho <- 2
# 
# df <- hitchip_pred_range 
# df <- df %>% mutate(sq = sigma*exp( -.5*(((x-y)/rho)^2)))
# df <- df %>% mutate(linear = 0.01*(x^2 + y^2))
# df <- df %>% mutate(sq.linear = sq - linear)
# 
# 
# df %>% ggplot(aes(x = x, y = y, fill = sq)) + geom_tile() + 
#   scale_fill_gradient2(low = "blue", high = "red")

## SDE stationary density ******************** ####
SDE_non_stationary_density <- function(X, drift, diffusion) {
  
  if((length(X) != length(drift)) |
     (length(X) != length(diffusion))) stop("invalid input lengths")
  
  dX <- X[2] - X[1]
  
  # Scale measure
  scale_measure <- exp(-2*cumsum(drift/diffusion^2)*dX)
  
  # Speed measure
  speed_measure <- 1/(scale_measure*diffusion^2)
  
  speed_measure <- speed_measure/(sum(speed_measure)*dX)
  
  
  return(data.frame(x = X, y = speed_measure))
  
}

