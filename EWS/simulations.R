## Simulate non-stationary AR(p) series ************ ####

simulate_ar <- function(N, p = 1, coefs, c, sigma) {
  
  if(length(coefs) != p) {
    stop("Invalid parameters")
  }
  
  x <- 1:N
  y <- rep(NA, N)
  
  # Random intial value
  y[1:p] <- rnorm(p, 0, 1)
  
  for(i in (p + 1):N) {
    
    y[i] <- next_ar_value(y[(i-1):(i-p)], p, coefs, c, sigma)
    
  }
    
  return(y)
    
}

next_ar_value <- function(prev, p = 1, coefs, c, sigma) {
  
  if(length(prev) != p) {
    stop("Too few previous observations")
  }
  
  next_y <- c + sum(coefs[1:p]*prev[p:1]) + rnorm(1, 0, sigma)
  
  return(next_y)
}

# AR(1): constant c, coefs[1] linearly from 0 to 1
N <- 1000
# coefs <- (0:(N-1))/N
coefs <- data.frame(ph1 = seq(from = .25, to = .75, length.out = 1000), 
                    ph2 = seq(from = 0, to = .25, length.out = 1000))
sigma <- seq(from = .25, to = 1, length.out = N)

simulate_non_stat_ar <- function(N, p = 1, coefs, c, sigma) {
  # 
  # if(length(coefs) != p) {
  #   stop("Invalid parameters")
  # }
  # 
  
  if(p == 1) {
    coefs <- coefs %>% as.data.frame()
  }
  
  x <- 1:N
  y <- rep(NA, N)
  
  # Random intial value
  y[1:p] <- rnorm(p, 0, 1)
  
  
  for(i in (p + 1):N) {
    
    y[i] <- next_ar_value(y[(i-1):(i-p)], p, coefs[i,], c, sigma)
    
  }
  
  return(y)
  
}

simulate_non_stat_ar(N, p = 2, coefs, 0, .1) %>% plot(type = "l")



# Alla olevan juuret z1, z2 pitää olla yksikköympyrän ulkopuolella, jotta stationaarinen
# 1 - phi1*z - phi2*z^2
# 
# 
# -phi2*(z - z1)*(z - z2) = 
#   -phi2*z^2 + phi2*(z1 + z2)*z - phi2*z1*z2
# 
# 
# -phi2*z1*z2 == 1
# phi2*(z1 + z2) == -phi1
# 
# phi1 == (z1 + z2)/(z1*z2)
# phi2 == -1/(z1*z2)


# Jos pidetään z2 paikallaan = .25 ja kasvatetaan z1 .5 --> 1
z2 <- rep(10, N)
z1 <-  seq(from = 5, to = 2, length.out = N)

phi1 <-  (z1 + z2)/(z1*z2)
phi2 <- -1/(z1*z2)

coefs <- data.frame(phi1, phi2)
simulate_non_stat_ar(N, p = 2, coefs, 0, .1) %>% plot(type = "l")




# data.frame(phi1, phi2) %>% 
#   ggplot(aes(x = phi1, y = phi2)) +
#   geom_path()

## Vostok ice core ********************* ####
vostok <- read.table(file = "EWS/data/vostok_ice_core.txt",
                     sep = "\t", header = TRUE)



vostok <- vostok %>% 
  set_colnames(c("Depth", "Age", "Deut", "Delta"))


# 
# vostok %>% 
#   mutate(Age = -Age) %>%
#   ggplot(aes(x = Age, y = Delta)) + 
#   geom_point()


# Add eras
vostok <- vostok %>% 
  mutate(Era = NA) %>% 
  mutate(Era = ifelse(Age < 58800 & Age > 12000, 1, Era)) %>% 
  mutate(Era = ifelse(Age < 151000 & Age > 128000, 2, Era)) %>% 
  mutate(Era = ifelse(Age < 270000 & Age > 238000, 3, Era)) %>% 
  mutate(Era = ifelse(Age < 385300 & Age > 324600, 4, Era))
  

# Add transition phases
vostok <- vostok %>% 
  mutate(Phase = ifelse(!is.na(Era), "pre", NA)) %>% 
  mutate(Phase = ifelse(Age > 12000 & Age < 17000, "post", Phase)) %>% 
  mutate(Phase = ifelse(Age > 128000 & Age < 135000, "post", Phase)) %>% 
  mutate(Phase = ifelse(Age > 238000 & Age < 242000, "post", Phase)) %>% 
  mutate(Phase = ifelse(Age > 324600 & Age < 334100, "post", Phase))
  
# Time scale to kilo years & Flip axis
vostok <- vostok %>% 
  mutate(Age = -Age/1000)
  
# Plot
vostok %>% 
  filter(!is.na(Era)) %>% 
  ggplot(aes(x = Age, y = Delta)) + 
  geom_line(aes(color = Phase)) +
  facet_wrap(~Era, scales = "free") +
  scale_color_tron()



interpolated_vostok <- lapply(1:4, function(i) {
  
  df <- vostok %>% 
    filter(Era == i, Phase == "pre")
    
  # min_diff <- (df$Age[2:dim(df)[1]] - df$Age[1:(dim(df)[1] - 1)]) %>% abs %>% mean
  
  df <- df %>% 
    mutate(Age = Age - min(Age)) %>% 
    mutate(Age = (Age/max(Age))*nrow(df)) %>% 
    arrange(Age)
  
  interpolated <- lapply(1:(nrow(df) - 1), function(x) {
    
    prev <- last(which(df$Age < x))
    cons <- first(which(df$Age > x))
    
    sub_df <- df[prev:cons, ]
    
    if(nrow(sub_df) == 3) {
      return(sub_df[2, "Delta"])
    } else {
      
      x1 <- sub_df[1, "Age"]
      x2 <- sub_df[2, "Age"]
      y1 <- sub_df[1, "Delta"]
      y2 <- sub_df[2, "Delta"]
      
      interpolated_value <- x*(y2 - y1)/(x2 - x1) + (x2*y1 - x1*y2)/(x2 - x1)
      
      return(interpolated_value)
    }
    
  }) %>% unlist
  
  interpolated <- c(df$Delta[1], interpolated)
  
  data.frame(Delta = interpolated, 
             Age = 1:length(interpolated), 
             Era = i, 
             Phase = "pre") %>% 
    return()
  
  
  
}) %>% do.call(rbind, .) %>% 
  mutate(interpolation = TRUE)

vostok_edit <-lapply(1:4, function(i) {
  
  df <- vostok %>% 
    filter(Era == i, Phase == "pre")
  
 
  df <- df %>% 
    mutate(Age = Age - min(Age)) %>% 
    mutate(Age = (Age/max(Age))*nrow(df)) %>% 
    select(-c(Depth, Deut))
  
  return(df)
  
}) %>% do.call(rbind, .) %>% 
  mutate(interpolation = FALSE)
  

vostok_pre <- rbind(interpolated_vostok, vostok_edit)

saveRDS(vostok_pre, file = "EWS/data/vostok_pre.RDS")

ggplot() + 
  geom_point(data = interpolated_vostok, aes(x = Age, y = Delta), color = "red") +
  geom_point(data = vostok_edit, aes(x = Age, y = Delta)) +
  facet_wrap(~Era, scales = "free")
