# Lambda per sigma

set.seed(1)
times <- seq(from = 0, to = 250, by = 0.1)
N_series <- 100

ews_set <- lapply(1:N_series, function(j) {
  
  r <- 1
  K <- 10
  cs <- 1
  h <- 1
  sigma <- 0.07
  
  
  y <- rep(NA, length(times))
  y[1] <- 8
  
  cs <- 1 + (2.7 - 1)*(times/max(times))
  
  
  for(i in 2:length(times)) {
    
    
    # c <- cs[i]

    dt <- times[(i-1):i]
    
    y[i] <- ews_generator(y[i-1], dt, c = cs[i], milstein = T)[2]

    
  }
  
  
  return(log(y + 1))

  
}) %>%
  do.call(cbind, .) %>%
  cbind(x = times) %>%
  as.data.frame()



# ews_set %>% as.data.frame() %>% melt(id.vars = "x") %>% ggplot(aes(x = x, y = value, color = variable)) + geom_line()



ews_set <- ews_set[1:(which((cs < 2.604)) %>% max()), ]

ews_set <- ews_set %>% thin


# Get EWS with generic_ews
ews_set_stats <- lapply(1:N_series, function(i) {
  
  stats <- generic_ews(ews_set[, i], winsize = 50) 
  
  stats <- stats %>% 
    mutate(resilience = ar1/(sd^2))
  
  
  dev.off()
  
  
  stats %>% 
    select(timeindex, ar1, sd, resilience) %>% 
    mutate(series = i)%>%
    as.data.frame()
  
})  %>% 
  do.call(rbind, .)



# Standadize
ews_set_stats <- ews_set_stats %>% 
  mutate(ar1 = (ar1 - mean(ar1))/sd(ar1), 
         sd = (sd - mean(sd))/sd(sd), 
         resilience = (resilience - mean(resilience))/sd(resilience))


# Compute difference of AR1, Var and AR1/Var in first and last windows
ews_set_stats_delta <- lapply(1:N_series, function(i) {
  
  # Select serires
  df <- ews_set_stats %>%
    filter(series == i)
  
  
  # First and last time points
  minimum <- which(min(df$timeindex) == df$timeindex)
  maximum <- which(max(df$timeindex) == df$timeindex)
  
  
  df <- (df[minimum, c("ar1", "sd", "resilience")] - df[maximum, c("ar1", "sd", "resilience")]) 
  
  df %>% mutate(series = i)
  
}) %>% do.call(rbind, .)


# Plot
ews_set_stats_delta %>% melt(id.vars = "series") %>%
  ggplot(aes(x = variable, y = abs(value))) +
  geom_point() + 
  geom_line(aes(group = series)) + 
  geom_boxplot()


# Proportion of series in which change of resilience is larger than that of ar1 and sd
ews_set_stats_delta %>%
  group_by(series) %>% 
  summarise(resilience > abs(ar1) & resilience > abs(sd)) %>% 
  set_colnames(c("series", "test")) %>% 
  ungroup %>% 
  pull(test) %>% mean

