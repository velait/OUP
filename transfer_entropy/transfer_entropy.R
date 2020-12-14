pseq <- stool_pseq_genus_clr %>% subset_samples(host == "M3")

df <- t(abundances(pseq)) %>% 
  cbind(time = meta(pseq)$days_since_experiment_start)

df <- df %>% 
  as.data.frame() %>% 
  arrange(time) 


all_pairs <- all_unique_pairs(taxa_names(pseq))

TE_pairs <- lapply(1:nrow(all_pairs), function(i) {
  
    print(i)
    TE <- transfer_entropy(df[, all_pairs[i, 1]], df[, all_pairs[i, 2]])
    
    coefs <- TE$coef
    
    res <- data.frame(x = all_pairs[i, c(1:2)], y = all_pairs[i, c(2:1)], 
                      coefs)
    
    
    res <- res %>% set_rownames(NULL)
    
    return(res)
}) %>% do.call(rbind, .)




all_unique_pairs <- function(x) {
  
  all_pairs <- lapply(1:(length(x) - 1), function(i) {
    lapply((i+1):length(x), function(j) {
      
      c(x[i], x[j])
      
    }) %>% do.call(rbind, .)
    
  }) %>% do.call(rbind, .)
  
  
  return(all_pairs)
}



TE_pairs %>% 
  mutate(significant = ifelse(p.value <= 0.05, "yes", "no")) %>% 
  ggplot(aes(x = x, y = y, fill = significant)) + 
  geom_tile()


TE_pairs %>% 
  mutate(significant = ifelse(p.value <= 0.05, "yes", "no")) %>% 
  filter(significant == "yes") %>% 
  ggplot(aes(x = x, y = y, fill = te)) + 
  geom_tile() + 
  scale_fill_gradientn(colours = c("white", "black", "red", "blue"))



xx <- TE_pairs %>% 
  mutate(significant = ifelse(p.value <= 0.05, "yes", "no")) %>% 
  filter(significant == "yes") %>% 
  arrange(rev(te)) %>% 
  select(x, y)


i <- 12
df[, c(xx[i+11, "x"] %>% 
         as.character(),
       xx[i, "y"] %>% 
         as.character(),
       "time")] %>% 
  as.data.frame() %>% 
  melt(id.vars = "time") %>% 
  ggplot(aes(x = time, y = value, color = variable)) +
  geom_line()
