# Sensitivity analysis for leading indicators in microbiome like data


## Window/Bandwidth analyzer  *************************** ####
sensitivity_analyzer <- function(data) {
  
  mydata <- data
  
  N <- ncol(mydata) - 1
  
  # windows
  res <- lapply(seq(from = 25, to = 75, by = 5), function(win) {
    
    print(paste("win", win))
    
    # bandwidth
    bw_res <- lapply(seq(from = 5, to = 30, by = 5), function(bw) {
      
      # data
      set_tau <- lapply(1:N, function(i) {
        
        # Get ews
        ews <- generic_ews(mydata[, i],
                           winsize = win, 
                           detrending = "gaussian",
                           bandwidth = bw)
        
        # Close plot window
        dev.off()
        
        # Kendall's tau
        tau <- cor(ews$ar1, ews$timeindex, method = "kendall")
        
      }) %>% unlist
      
      return(data.frame(tau = set_tau, window = win, bandwidth = bw, series = 1:N))
      
    }) %>%
      do.call(rbind, .)
    
    
  }) %>%
    do.call(rbind, .)
  
  return(res)
  
}


set_res_dt_1 <- sensitivity_analyzer(ews_set)

set_res_dt_0.1 <- sensitivity_analyzer(ews_set)


set_res %>% 
  group_by(window, bandwidth) %>%
  summarize(mean(tau)) %>%
  set_colnames(c("window", "bandwidth", "tau")) %>%
  ggplot(aes(x = window, y = bandwidth, fill = tau)) +
  geom_tile()  +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue")
  



## Delta length scale analyzer  ************************* ####


length_scale_analyzer <- function(N) {
  
  
  # Initial length scale
  res <- lapply(seq(from = 1, to = 30, by = 1), function(rho1) {
    
    # Length scale change
    lapply(seq(from = 1, to = 30, by = 1), function(delta_rho) {
      
      
      # Generate data
      x <- seq(from = 0, to = 100, length.out = 101)
      sigma <- seq(from = .5, to = .75, length.out = length(x))
      rho <- seq(from = rho1, to = rho1 + delta_rho, length.out = length(x))
      epsilon <- .01
      
      mydata <- oup_input_set(N, x, rho, sigma, epsilon)
      
      
      # data
      set_tau <- lapply(1:N, function(i) {
        
        # Get ews
        ews <- generic_ews(mydata[, i])
        
        # Close plot window
        dev.off()
        
        # Kendall's tau
        tau <- cor(ews$ar1, ews$timeindex, method = "kendall")
        
      }) %>% unlist
      
      data.frame(tau = set_tau, series = 1:N, rho1 = rho1, delta_rho = delta_rho)
      
    }) %>% do.call(rbind, .)
    
  }) %>% do.call(rbind, .)
  
  
  return(res)
}

length_scale_analyzer2 <- function(N,
                                   x = seq(from = 0, to = 100, length.out = 101), 
                                   sigma = seq(from = .25, to = .5, length.out = 101), 
                                   rho = seq(from = 1, to = 20, by = .25), 
                                   delta_rho = seq(from = 1, to = 20, by = 0.25)) {
  
  
  # Initial length scale
  res <- lapply(rho, function(r) {
    
    # Length scale change
    lapply(delta_rho, function(delta_r) {
      
      
      # Generate data
      # x <- seq(from = 0, to = 100, length.out = 101)
      # sigma <- seq(from = .25, to = .5, length.out = length(x))
      rho <- seq(from = r, to = r + delta_r, length.out = length(x))
      autoc <- exp(-1/rho)
      epsilon <- .01

      # mydata <- oup_input_set(N, x, rho, sigma, epsilon)
      mydata <- square_exp_input_set(N, x, rho, sigma, epsilon)
      
      
      # data
      set_tau <- lapply(1:N, function(i) {
        
        # Get ews
        ews <- generic_ews(mydata[, i])
        
        # Close plot window
        dev.off()
        
        # Kendall's tau
        tau <- cor(ews$ar1, ews$timeindex, method = "kendall")
        
      }) %>% unlist
      
      data.frame(tau = set_tau, series = 1:N,
                 rho1 = rho[1], delta_rho = rho[length(rho)] - rho[1])
      
    }) %>% do.call(rbind, .)
    
  }) %>% do.call(rbind, .)
  
  
  # Add sigma range as character
  
  res <- res %>% 
    mutate(sigma = paste0(sigma[1], " - ", sigma[length(sigma)]))
  
  
  return(res)
}




sensi_rho1 <- length_scale_analyzer2(N = 1,
                                     x = 1:201,
                                    sigma = seq(from = .25, to = .5, length.out = 201),
                                    rho = seq(from = 1, to = 10, by = .5),
                                    delta_rho = seq(from = 1, to = 10, by = .5))

sensi_rho2 <- length_scale_analyzer2(N = 1,
                                     x = 1:201,
                                    sigma = seq(from = .25, to = .5, length.out = 201),
                                    rho = seq(from = 1, to = 10, by = .5),
                                    delta_rho = seq(from = 10, to = 20, by = .5))

sensi_rho3 <- length_scale_analyzer2(N = 1,
                                     x = 1:201,
                                    sigma = seq(from = .25, to = .5, length.out = 201),
                                    rho = seq(from = 10, to = 20, by = .5),
                                    delta_rho = seq(from = 1, to = 10, by = .5))

sensi_rho4 <- length_scale_analyzer2(N = 1,
                                     x = 1:201,
                                    sigma = seq(from = .25, to = .5, length.out = 201),
                                    rho = seq(from = 10, to = 20, by = .5),
                                    delta_rho = seq(from = 10, to = 20, by = .5))


sensi_rho <- rbind(sensi_rho1, 
                   sensi_rho2, 
                   sensi_rho3, 
                   sensi_rho4)


 p_high <- sensi_rho %>% 
  # group_by(delta_rho, rho1, sigma) %>% 
  # summarize(mean_tau = mean(tau)) %>%
  ggplot() +
  geom_tile(aes(x = rho1, y = delta_rho, fill = tau)) +
  # geom_raster(aes(x = rho1, y = delta_rho, fill = tau),
  #             interpolate = TRUE, ) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  facet_wrap(.~sigma, labeller = labeller(sigma)) + 
  labs(subtitle = "Average tau over 20 series in each tile")


sensi_rho1 %>% filter(rho1 == 1, delta_rho == 1)
  ggplot() +
  geom_tile(aes(x = rho1, y = delta_rho, fill = tau)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) 



  plot_grid(p_low + labs(subtitle = "100 time points"), p_high + labs(subtitle = "200 time points"))
  
  
xx %>% 
  group_by(ac1, ac2) %>% 
  summarize(mean_tau = mean(tau)) %>%
  ggplot() +
  geom_point(aes(x = ac1, y = ac2, color = mean_tau)) +
  geom_abline(slope = 1, intercept = 0) + 
  # geom_raster(aes(x = rho1, y = delta_rho, fill = mean_tau),
  #             interpolate = TRUE) + 
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-.5, .5)) +
  labs(subtitle = "Average tau over 20 series in each tile")
