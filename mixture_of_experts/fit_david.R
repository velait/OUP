# Fit david et al to mixture of experts
library(microbiome)

load("/Users/villelaitinen/Desktop/PhD/early_warning_signals/data-David2014/David_phyloseq.Rdata")

seqA.clr <- seqA %>% aggregate_taxa("Genus") %>%  microbiome::transform("compositional") %>% microbiome::transform("log10")
seqA.clr.genus <- seqA.clr %>% abundances() %>% t %>% as.data.frame() %>% mutate(time=meta(seqA)$time)

seqB.clr <- seqB %>% aggregate_taxa("Genus") %>%  microbiome::transform("compositional") %>%  microbiome::transform("log10")
seqB.clr.genus <- seqB.clr %>% abundances() %>% t %>% as.data.frame() %>% mutate(time=meta(seqB)$time)



bimodality.score <- bimodality(seqA.clr, method = "potential_analysis",
                               bs.iter = 10, peak.threshold = 10,
                               min.density = 10)

bimodality.score <- bimodality.score %>%
  as.data.frame() %>% 
  rownames_to_column() %>% 
  set_colnames(c("genus", "coef")) %>% 
  arrange(desc(coef))

genus <- bimodality.score$genus[8]

# plot time series
example_series <- seqB.clr.genus[, c(genus, "time")] %>% 
  set_colnames(c("genus", "time"))

example_series %>%   
ggplot(aes(x = time, y = log(genus)))+
geom_point(color = "black") +
geom_line()


# Histogram
# seqB.clr.genus[, c(genus)] %>% hist(breaks = 100)

## Stan

# data
david_data <- list(y = seqB.clr.genus[, c(genus)], 
     x = seqB.clr.genus[, c("time")],
     N = length(seqB.clr.genus[, c("time")]))


# fit
david_samples <- sampling(fit_oup_mixture,
                                david_data,
                                chains = 1,
                                iter = 1000, 
                                init = -10)



# Plot the inferred latent function
david_fit_weights <- summary(david_samples)$summary[grep("weight", rownames(summary(david_samples)$summary)), c("mean","25%", "50%", "75%")]

cbind(x = david_data$x, david_fit_weights) %>% 
  set_colnames(c("x", "mean", "lower", "mode", "upper")) %>% 
  as.data.frame() %>% 
  ggplot() +
  geom_line(aes(x = x, y = mode)) +
  geom_ribbon(aes(x = x, ymin = lower, ymax = upper), fill = "grey", alpha = .25) +
  theme_bw() +
  geom_line(data = example_series, aes(x = time, y = genus))
  


