load(file = "/Users/villelaitinen/Desktop/PhD/early_warning_signals/data-Ravel2012/data.clr.rda")

# 32 subjects
# 95 time points
# 217 genera

data.cube %>% class


ravel_df <- lapply(1:dim(data.cube.clr)[1], function(i) data.cube.clr[i, , ])

ravel_df_combined <- ravel_df %>% do.call(rbind, .)

pca <- prcomp(ravel_df_combined)

pca_clr_time <- cbind(pca$x[, 1:2],
                      time =  as.numeric(ravel_df_combined %>% rownames),
                      subject = rep(1:32, each = 95)) %>% 
  as.data.frame()

rownames(pca_clr_time) <- NULL

pca_clr_time %>%
  filter(subject == 2) %>% 
  as.data.frame() %>% 
  ggplot(aes(x = PC1, y = PC2)) + 
  geom_path() +
  facet_wrap(~subject)


write.csv(pca_clr_time %>% filter(subject == 2),
          file = "NPDE/ravel_df.csv", row.names = FALSE)

write.csv(pca_clr_time %>% filter(subject == 3),
          file = "NPDE/ravel_df2.csv", row.names = FALSE)



pca_clr_time %>%
  filter(subject == 2) %>% 
  as.data.frame() %>% 
  ggplot(aes(x = time, y = PC1)) + 
  geom_point()
  