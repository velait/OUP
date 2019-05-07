C <- 27*alpha^2 - 4*beta^3

C_df <- lapply(seq(-2, 2, length.out = 100), function(a) {
  
  sapply(seq(-2, 2, length.out = 100), function(b) {
    
    c(a = a, b = b, C = 27*a^2 - 4*b^3)
    
  } ) %>% t
  
}) %>%
  do.call(rbind, .) %>% 
  as.data.frame()



C_df %>% 
  mutate(sign = ifelse(C < 0 , 0, 1)) %>% 
  ggplot() +
  geom_tile(aes(x = a, y = b, fill = C)) +
  geom_contour(aes(x = a, y = b, z = sign), color = "black", binwidth = 1) + 
  scale_fill_gradient2(low = "blue", high = "red", mid = "white") +
  labs(title = "Cardano's discriminant",
       x = expression(alpha), y = expression(beta),
       fill = "")
  

