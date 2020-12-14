library(MicrobeDS)
library(microbiome)
library(magrittr)
library(tidyverse)
library(reshape2)
library(earlywarnings)
library(changepoint.np)

## Data ********************************** ####

data("MovingPictures")

# Filter relevnat meta data
sample_data(MovingPictures) <- meta(MovingPictures)[, c("days_since_experiment_start", "host", "sample_type", "common_sample_site", "days_since_epoch")]


# Separate phyloseqs from different sample types
stool_pseq <- subset_samples(MovingPictures, common_sample_site == "feces")
hand_L_pseq <- subset_samples(MovingPictures, common_sample_site == "L_palm")
hand_R_pseq <- subset_samples(MovingPictures, common_sample_site == "R_palm")
saliva_pseq <- subset_samples(MovingPictures, common_sample_site == "Tongue")


# Get rare OTUs in each sample type
super_rare_otus <- list(hand_L = names(which((abundances(hand_L_pseq) %>% rowSums) > 10)),
                        hand_R = names(which((abundances(hand_R_pseq) %>% rowSums) > 10)),
                        stool = names(which((abundances(stool_pseq) %>% rowSums) > 1)), 
                        saliva = names(which((abundances(saliva_pseq) %>% rowSums) > 1)))


# Filter out rare OTUs
stool_pseq@otu_table <- stool_pseq@otu_table[super_rare_otus[["stool"]], ]
hand_L_pseq@otu_table <- hand_L_pseq@otu_table[super_rare_otus[["hand_L"]], ]
hand_R_pseq@otu_table <- hand_R_pseq@otu_table[super_rare_otus[["hand_R"]], ]
saliva_pseq@otu_table <- saliva_pseq@otu_table[super_rare_otus[["saliva"]], ]


# Aggregate to genus level
stool_pseq_genus <- aggregate_taxa(stool_pseq, level = "Genus")
hand_L_pseq_genus <- aggregate_taxa(hand_L_pseq, level = "Genus")
hand_R_pseq_genus <- aggregate_taxa(hand_R_pseq, level = "Genus")
saliva_pseq_genus <- aggregate_taxa(saliva_pseq, level = "Genus")

# saveRDS(stool_pseq_genus, file = "data/moving_pictures/stool_pseq_genus.Rdata")
# saveRDS(hand_L_pseq_genus, file = "data/moving_pictures/hand_L_pseq_genus.Rdata")
# saveRDS(hand_R_pseq_genus, file = "data/moving_pictures/hand_R_pseq_genus.Rdata")
# saveRDS(saliva_pseq_genus, file = "data/moving_pictures/saliva_pseq_genus.Rdata")

stool_pseq_genus <- readRDS(file = "data/moving_pictures/stool_pseq_genus.Rdata")
hand_L_pseq_genus <- readRDS(file = "data/moving_pictures/hand_L_pseq_genus.Rdata")
hand_R_pseq_genus <- readRDS(file = "data/moving_pictures/hand_R_pseq_genus.Rdata")
saliva_pseq_genus <- readRDS(file = "data/moving_pictures/saliva_pseq_genus.Rdata")

# Compositional transformation
saliva_pseq_genus_compositional <- transform(saliva_pseq_genus, "compositional")
hand_L_pseq_genus_compositional <- transform(hand_L_pseq_genus, "compositional")
hand_R_pseq_genus_compositional <- transform(hand_R_pseq_genus, "compositional")
stool_pseq_genus_compositional <- transform(stool_pseq_genus, "compositional")


saliva_pseq_genus_core <- core(saliva_pseq_genus_compositional, detection = .1/100, prevalence = 10/100)
hand_L_pseq_genus_core <- core(hand_L_pseq_genus_compositional, detection = .1/100, prevalence = 10/100)
hand_R_pseq_genus_core <- core(hand_R_pseq_genus_compositional, detection = .1/100, prevalence = 10/100)
stool_pseq_genus_core <- core(stool_pseq_genus_compositional, detection = .1/100, prevalence = 10/100)


# # CLR transformation
saliva_pseq_genus_clr <- transform(saliva_pseq_genus_core, "clr")
hand_L_pseq_genus_clr <- transform(hand_L_pseq_genus_core, "clr")
hand_R_pseq_genus_clr <- transform(hand_R_pseq_genus_core, "clr")
stool_pseq_genus_clr <- transform(stool_pseq_genus_core, "clr")

# Log10 transformation
# saliva_pseq_genus_log <- transform(saliva_pseq_genus_core, "log10")
# sebum_pseq_genus_log <- transform(sebum_pseq_genus_core, "log10")
# stool_pseq_genus_log <- transform(stool_pseq_genus_core, "log10")



## Hclust plot *************************** ####

hclust_tile_plot <- function(pseq, na_rm = FALSE) {
  
  mypseq <- pseq
  
  
  mydf <- cbind(t(abundances(mypseq)) %>% 
                  as.data.frame(), 
                meta(mypseq) %>% 
                  select(days_since_experiment_start))
  
  mydf <- mydf %>% 
    arrange(days_since_experiment_start)
  
  df <- mydf %>% select(-days_since_experiment_start)
  
  
  df_clust <- hclust(dist(t(df), method = "euclidean"), method = "ward.D")
  
  
  
  df_data <- df[, df_clust$order]
  
  p <- df_data %>% 
    cbind(time = mydf$days_since_experiment_start, .) %>%
    # cbind(time = 1:nrow(df_data)) %>% 
    melt(id.vars = "time") %>% 
    ggplot(aes(x = time, y = variable, fill = value) ) +
    geom_tile() + 
    scale_fill_gradientn(colors = c("white", "black", "red"), na.value = "grey")
  
  
  return(p)
}
hclust_tile_difference_plot <- function(pseq, square = FALSE, na_rm = FALSE) {
  
  mypseq <- pseq
  
  
  mydf <- cbind(t(abundances(mypseq)) %>% 
                  as.data.frame(), 
                meta(mypseq) %>% 
                  select(days_since_experiment_start))
  
  mydf <- mydf %>% 
    arrange(days_since_experiment_start)
  
  df <- mydf %>% select(-days_since_experiment_start)
  
  
  df_clust <- hclust(dist(t(df), method = "euclidean"), method = "ward.D")
  
  
  
  df_data <- cbind(time = mydf$days_since_experiment_start, df[, df_clust$order])
  
  diff_data <- (df_data[2:nrow(df_data), ] - df_data[1:(nrow(df_data) - 1), ]) %>% 
    mutate(time = time + df_data[1:(nrow(df_data) - 1), "time"])
  
  if(square) diff_data[, -1] <- diff_data[, -1]^2
  
  p <- diff_data %>% 
    # cbind(time = 1:nrow(df_data)) %>% 
    melt(id.vars = "time") %>% 
    ggplot(aes(x = time, y = variable, fill = value) ) +
    geom_tile() + 
    scale_fill_gradientn(colors = c("white", "black", "blue", "red"), na.value = "grey")
  
  
  return(p)
}
plot_otu_time_series <- function(pseq, otu, xlim = "all") {
  
  mypseq <- pseq
  
  
  mydf <- data.frame(t(abundances(mypseq))[, otu] %>% 
                       as.data.frame(), 
                     meta(mypseq) %>% 
                       select(days_since_experiment_start)) %>% 
    set_colnames(c(otu, "time"))
  
  if(class(xlim) == "numeric") {
    mydf <- mydf %>% 
      filter(time >= xlim[1] & time < xlim[2])
  }
  
  
  p <- mydf %>%
    melt(id.vars = "time") %>% 
    ggplot(aes(x = time, y = value)) + 
    geom_line() + 
    geom_point(size = .5) +
    facet_wrap(~variable, scales = "free")
  
  
  
  return(p)
}
pseq_generic_ews <- function(pseq, otu, xlim = "all", detrending = "gaussian", window = 50) {
  
  mypseq <- pseq
  
  
  mydf <- data.frame(t(abundances(mypseq))[, otu] %>% 
                       as.data.frame(), 
                     meta(mypseq) %>% 
                       select(days_since_experiment_start)) %>% 
    set_colnames(c(otu, "time"))
  
  if(class(xlim) == "numeric") {
    mydf <- mydf %>% 
      filter(time >= xlim[1] & time < xlim[2])
  }
  
  
  generic_ews(mydf[, otu], detrending = detrending, winsize = window)
  
}


pseq <- hand_L_pseq_genus_clr %>%
  subset_samples(host == "M3")

p <- hclust_tile_plot(pseq)
p


q <- hclust_tile_difference_plot(pseq, square = TRUE)
q

taxon <- c("4430826")

plot_otu_time_series(pseq,
                     otu = taxon, 
                     xlim = "all") +
  geom_vline(xintercept = c(280))




pseq_generic_ews(pseq,
                 otu = otu, 
                 xlim =  c(0, 280), 
                 window = 10)

## Change points ************************* ####
pseq_change_point_detector <- function(pseq, otu, min_seq_length = 10, rank = FALSE) {
  
  mypseq <- pseq
  
  
  
  mydf <- data.frame(t(abundances(mypseq))[, otu] %>% 
                       as.data.frame(), 
                     meta(mypseq) %>% 
                       select(days_since_experiment_start)) %>% 
    set_colnames(c(otu, "time")) %>% 
    arrange(time)
  
  
  cp_res <- lapply(1:length(otu), function(i) {
    
    x <- mydf[, c(otu[i], "time")]
    
    cp <- cpt.np(x[, 1], minseglen = min_seq_length)@cpts
    
    if(isTRUE(rank)) {
      change_point <- cp[1:(length(cp) - 1)]
    } else {
      change_point <- x[cp[1:(length(cp) - 1)], "time"]
    }
    
    
    return(data.frame(change_point = change_point, otu = otu[i]))
    
  }) %>% do.call(rbind, .)
  
  
  
  return(cp_res)
  
}

pseq_cp <- pseq_change_point_detector(pseq,
                                      otu = taxa_names(pseq), min_seq_length = 20, rank = FALSE)



p +
  geom_tile(data = pseq_cp,
            aes(x = change_point, y = otu),
            fill = "green", na.rm = T) + 
  labs(title = "Moving pictures; M3, left hand; green = change point")


taxon <- taxa_names(pseq)[1]

plot_otu_time_series(pseq,
                     otu = taxon,
                     xlim = "all") +
  geom_vline(data = (pseq_cp %>%
                       filter(otu %in% taxon) %>%
                       mutate(variable = otu)), aes(xintercept = change_point), linetype = "dashed", color = "red")

## EWS *********************************** ####

# Segment time series according to time points starters, otu is added as a column in the results matrix
segment_ts <- function(x, time, starters, otu) {
  
  # segment starting times. Add first and last time index to aid segmentation. 
  starters <- c(time[1], starters, time[length(time)]) %>% unname
  starters <- starters[!duplicated(starters)]
  
  # starting indices
  starter_ind <- which(time %in%  starters)
  
  # Segmentation
  seg_res <- lapply(1:(length(starter_ind) - 1), function(i) {
    
    seg <- x[starter_ind[i]:starter_ind[i + 1]]
    seg_time <- time[starter_ind[i]:starter_ind[i + 1]]
    seg_length <- time[starter_ind[i + 1]]- time[starter_ind[i]]
    seg_n <- length(seg_time)
    seg_start <- starters[i]
    seg_end <- max(time[time < starters[i + 1]])
    
    return(data.frame(seg, seg_time, i, otu, seg_length, seg_n, seg_start, seg_end) %>%
             set_colnames(c("abundance", "time", "segment", "otu",
                            "segment_length", "segment_n", "segment_start", "segment_end")))
  }) %>% do.call(rbind, .)
  
  
  return(seg_res) 
}

# Generic ews results according to earlywarnings package 
get_pseq_ews <- function(pseq,
                         otu = "all", 
                         min_seq_length = 20,
                         rank = FALSE) {
  
  
  mypseq <- pseq
  
  if(otu == "all") otu <- taxa_names(mypseq)
  
  # Get change points
  cp_df <- pseq_change_point_detector(mypseq,
                                      otu = otu,
                                      min_seq_length = min_seq_length,
                                      rank = rank)
  
  
  # Get time and abundances
  otu_df <- cbind(time = meta(mypseq)$days_since_experiment_start, 
                  t(abundances(mypseq))[, otu])
  
  # Arrange according to time column
  otu_df <- otu_df[order(otu_df[, "time"]), ]
  
  
  
  # Segment data
  seg_otu_df <- lapply(otu, function(i) {
    
    otu_cps <- cp_df %>%
      filter(otu == i) %>% 
      pull(change_point)
    
    
    otu_segments <- segment_ts(otu_df[, i],
                               otu_df[, "time"],
                               starters = otu_cps,
                               otu = i)
    
    return(otu_segments)
  }) %>%
    do.call(rbind, .)
  
  # EWS
  ews_res <- lapply(otu, function(i) {
    # print(i)
    mydf <- seg_otu_df %>% 
      filter(otu == i)
    
    otu_ews_res <- lapply(unique(mydf$segment), function(s) {
      
      seg_df <- mydf %>%
        filter(segment == s)
      
      segment_start <- seg_df$time[1]
      
      ews <-  seg_df %>% 
        pull(abundance) %>% 
        generic_ews(detrending = "gaussian") %>% 
        mutate(segment = s,
               otu = i,
               timeindex = as.numeric(timeindex) + segment_start, 
               segment_length = seg_df$segment_length[1], 
               segment_n = seg_df$segment_n[1], 
               segment_start = seg_df$segment_start[1], 
               segment_end = seg_df$segment_end[1])
      
      # Kendall's tau correlations
      tauski <- ews %>%
        melt(id.vars = c("timeindex", "otu", "segment",
                         "segment_length", "segment_n", "segment_start", "segment_end")) %>%
        mutate(value = as.numeric(value)) %>% 
        group_by(variable) %>% 
        summarize(tau = cor(timeindex, value, method = "kendall")) %>% 
        mutate(segment = s,
               otu = i, 
               segment_length = seg_df$segment_length[1], 
               segment_n = seg_df$segment_n[1], 
               segment_start = seg_df$segment_start[1], 
               segment_end = seg_df$segment_end[1])
        
      
      return(list(cor = tauski, par_df = ews))
      
    })
    
    # Remove pop-up plots
    for(d in dev.list()[-1]) dev.off(which = d)
    
    return(otu_ews_res)
    
  })
  
  # Combine tau and parameter trajectory lists
  tau_res <- do.call(rbind, lapply(ews_res, function(s) {
    lapply(s, function(sr) {
      
      sr[["cor"]]
      
    }) %>% do.call(rbind, .)
  } ))
  par_res <- do.call(rbind, lapply(ews_res, function(s) {
    lapply(s, function(sr) {
      
      sr[["par_df"]]
      
    }) %>% do.call(rbind, .)
  } ))
  
  
  return(list("tau" = tau_res, "parameters" = par_res))
  
} 


ews_res <- get_pseq_ews(pseq, 
                        otu = "all",
                        # otu = taxa_names(pseq)[1:2],
                        min_seq_length = 20, 
                        rank = FALSE)


# ews_res[["tau"]] %>% 
#   ggplot(aes(x = segment_n, y = tau)) + 
#   geom_point() +
#   geom_smooth(method = "loess") + 
#   facet_wrap(~variable)
  


tau_tile_plot <- function(pseq, ews_res) {
  
  mypseq <- pseq
  
  mydf <- cbind(t(abundances(mypseq)) %>% 
                  as.data.frame(), 
                meta(mypseq) %>% 
                  select(days_since_experiment_start))
  
  mydf <- mydf %>% 
    arrange(days_since_experiment_start)
  
  abundance_df <- mydf %>% select(-days_since_experiment_start)
  
  hclustering <- hclust(dist(t(abundance_df), method = "euclidean"), method = "ward.D")
  clustered_taxa <- colnames(abundance_df[, hclustering$order])
  df_data <- cbind(time = mydf$days_since_experiment_start, abundance_df[, hclustering$order])
  
  
  # Get tau correlations
  tau_res <- ews_res[["tau"]] %>% filter(variable %in% c("ar1", "sd"))
  
  
  # Color abundances based on ar1 and sd tau's
  tau_df <- lapply(c("ar1", "sd"), function(par) {
    
    par_tau_df <- lapply(taxa_names(mypseq), function(taxon) {
      
      otu_tau <- tau_res %>% filter(otu == taxon, variable == par)
      otu_data <- df_data[, c("time", taxon)]
      
      for(s in otu_tau$segment) {
        otu_tau_segment <- otu_tau %>% filter(segment == s)
        otu_data[otu_data$time >= otu_tau_segment$segment_start & 
                   otu_data$time <= otu_tau_segment$segment_end, taxon] <- otu_tau_segment$tau
      }
      
      return(otu_data[, taxon])
      
    }) %>% do.call(cbind, .) %>% set_colnames(taxa_names(mypseq))
    
    return(cbind(time = df_data$time, par_tau_df) %>% as.data.frame())
    
  }) %>% set_names(c("ar1", "sd"))
  
 
  plots <- lapply(names(tau_df), function(par) {
    df <- tau_df[[par]]
    
    p <- df[1:(nrow(df) - 1), c("time", clustered_taxa)] %>% 
      melt(id.vars = "time") %>% 
      ggplot(aes(x = time, y = variable, fill = value)) + 
      geom_tile() +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-1, 1)) +
      labs(x = "Time", y = "", subtitle = paste0("Time correlation, parameter ", par))
    return(p)
  }) %>% set_names(names(tau_df))
   
  
  return(plots)
}

parta_plot <- function(pseq, ews_res) {
  
  mypseq <- pseq
  
  mydf <- cbind(t(abundances(mypseq)) %>% 
                  as.data.frame(), 
                meta(mypseq) %>% 
                  select(days_since_experiment_start))
  
  mydf <- mydf %>% 
    arrange(days_since_experiment_start)
  
  abundance_df <- mydf %>% select(-days_since_experiment_start)
  
  hclustering <- hclust(dist(t(abundance_df), method = "euclidean"), method = "ward.D")
  clustered_taxa <- colnames(abundance_df[, hclustering$order])
  df_data <- cbind(time = mydf$days_since_experiment_start, abundance_df[, hclustering$order])
  
  
  # Get parameter trajectories
  par_trajectory <- ews_res[["parameters"]]
  

    
  plots <- lapply(c("ar1", "sd"), function(par) {
    
    df <- par_trajectory %>% select(timeindex, one_of(par), otu) %>% 
      set_colnames(c("Time", "Parameter", "otu"))
    
    p <- df %>% 
      ggplot(aes(x = Time, y = otu, fill = Parameter)) +
      geom_tile() + 
      scale_fill_gradientn(colours = c("red", "yellow", "green", "blue")) +
      labs(x = "Time", subtitle = paste0("Parameter trajectories, par = ", par))
    
    return(p)
  }) %>% set_names(c("ar1", "sd"))
  
  
  return(plots)
  
  
}

tau_tile_plot(pseq, ews_res)
parta_plot(pseq, ews_res)

## Potential analysis ******************** ####


# Potential analysis
bimodality_coefs <- lapply(list(saliva_pseq_genus_log,
                                sebum_pseq_genus_log, 
                                stool_pseq_genus_log), function(pseq) {
                                  
                                  bi_res <- list(M3 = bimodality(pseq %>% 
                                                                   subset_samples(host == "M3"),
                                                                 method = "potential_analysis"), 
                                                 F4 = bimodality(pseq %>% 
                                                                   subset_samples(host == "F4"),
                                                                 method = "potential_analysis"), 
                                                 Both = bimodality(pseq,
                                                                   method = "potential_analysis"))
                                  
                                  return(bi_res)
                                }) %>%
  set_names(c("saliva","sebum", "stool"))



abundances(saliva_pseq_genus_log)[which(bimodality_coefs[["saliva"]][["M3"]] > .6) %>% 
                                    names, ] %>% 
  t %>% 
  as.data.frame %>% 
  melt %>%
  ggplot(aes(x = value)) + 
  stat_density(geom = "line") + 
  facet_wrap(~variable, scales = "free")



# M3 bimodal plots

m3_df <- cbind(t(abundances((saliva_pseq_genus_log %>% 
                               subset_samples(host == "M3"))))%>% 
                 as.data.frame(), 
               meta(saliva_pseq_genus_log %>% 
                      subset_samples(host == "M3")) %>% 
                 select(days_since_experiment_start))

m3_df %>% 
  arrange(days_since_experiment_start) %>% 
  cbind(dummy = 1:nrow(m3_df)) %>% 
  select(-days_since_experiment_start) %>% 
  melt(id.vars = "dummy") %>%
  set_colnames(c("time", "variable", "value")) %>% 
  ggplot(aes(x = time, y = variable, fill = value)) +
  geom_tile() +
  scale_fill_gradientn(colors = c("blue", "green", "yellow", "red"), na.value = "grey")
# scale_fill_gradient2(midpoint = 0.025, low ="blue", mid = "white", high = "red")


# Potential examples
# M3 stool;
# GOOD: "4481131" @225, "356138" @225, 
# BAD: "2407149" @200, "4481427" @150, "163243" @120, "4454586" @300
# NON-STANDARD: "4306262", "215097", "4472174
# 
# 
# F4 stool; 
# Yes: "2407149" @60, 4465907@50
# Maybe: "3931537"
# L-V type dynamics (!): "1111115", "215097"
# 
# M3 Sebum;
# IFFY: "539293"



otu <- "4481131"

m3_df[, c("days_since_experiment_start", otu)] %>% 
  set_colnames(c("time", "abundance")) %>% 
  arrange(time) %>% 
  ggplot(aes(x = time, y = abundance)) +
  geom_point() +
  geom_line() +
  geom_vline(xintercept = 75, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 225, linetype = "dashed", color = "red")




m3_df[, c("days_since_experiment_start", otu)] %>% 
  set_colnames(c("time", "abundance")) %>% 
  filter(time >= 74 & time <= 225) %>%
  pull(abundance) %>% 
  # generic_ews()
  generic_ews(detrending = "gaussian", winsize = 50)




m3_df[, c("days_since_experiment_start", "1111115", "215097")] %>% 
  set_colnames(c("time", "otu1", "otu2")) %>% 
  filter(time >= 0 & time <= 117) %>% 
  melt(id.vars = "time") %>% 
  ggplot(aes(x = time, y = value, color = variable)) + 
  geom_line() + 
  geom_smooth(span = .3) +
  geom_point()











## PCA

df <- m3_df %>% arrange(days_since_experiment_start) %>% select(-days_since_experiment_start)


df_pca <- prcomp(df)


df_pca <- df_pca$x[, 2:3] %>% 
  as.data.frame() %>% 
  cbind(time = 1:nrow(df)) %>% 
  mutate(era = cut(time, breaks = 4))


df_pca %>% 
  melt(id.vars = c("time", "era")) %>% 
  # filter(era == unique(df_pca$era)[1]) %>% 
  ggplot() + 
  geom_line(aes(x = time, y = value, color = variable))

