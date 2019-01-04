# TODO #
# one series time series example: check submission?
# one series posteriors: example_plots.R
# potential well pictures: one and two wells
# hierarchical: how does things get better for a single one?
# none vs partial vs complete pooling?
# prior vs posterior
# two well sample time series

#### OUP graph ####


# lambda = 0.01, 0.1, 1
oup_example_sigma <- c(0.1, 0.25, 0.5)
oup_example_plot1 <- sapply(oup_example_sigma, function(x) generate_student_set(n_series=1, sigma = x, lambda = 0.5, mu = mu, intervals = 1:200, student_df=7, seed = 1)[["Y"]]) %>% as_tibble() %>% set_colnames(c("lambda = 0.01", "lambda = 0.1", "lambda = 1")) %>%  mutate(time=1:200) %>% melt(id.vars="time") %>% ggplot(aes(x=time, y=value, color=variable)) + geom_line() + scale_y_continuous(limits = c(-1, 1))  +  theme_bw() + labs(x="Time", y="", title="mu = 5, sigma = 0.1") + scale_color_manual(values=c('#377eb8', '#e41a1c', '#4daf4a')) + theme(legend.title=element_blank()) + ggtitle(TeX("$\\mu = 0;  \\lambda = 0.5$")) + scale_color_discrete(name = "$\\sigma$", labels = unname(c(TeX("$\\sigma = 0.1$"), TeX("$\\sigma = 0.25$"), TeX("$\\sigma = 0.5$"))))


# lambda = 0.01, 0.1, 1
oup_example_lambda <- c(0.1, 0.5, 1)
oup_example_plot2 <- sapply(oup_example_lambda, function(x) generate_student_set(n_series=1, sigma = 0.2, lambda = x, mu = mu, intervals = 1:200, student_df=7, seed = 1)[["Y"]]) %>% as_tibble() %>% set_colnames(c("lambda = 0.01", "lambda = 0.1", "lambda = 1")) %>%  mutate(time=1:200) %>% melt(id.vars="time") %>% ggplot(aes(x=time, y=value, color=variable)) + geom_line() + scale_y_continuous(limits = c(-1, 1)) +  theme_bw() + labs(x="Time", y="", title="mu = 5, sigma = 0.1") + scale_color_manual(values=c('#4daf4a','#377eb8', '#e41a1c')) + theme(legend.title=element_blank()) + ggtitle(TeX("$\\mu = 0;  \\sigma = 0.2$")) + scale_color_discrete(name = "$\\sigma$", labels = unname(c(TeX("$\\lambda = 0.1$"), TeX("$\\lambda = 0.5$"), TeX("$\\lambda = 1$"))))

#### One series inference ####

## Data
# sc_one_series_data <- generate_student_set(1, student_df = 7, mu = 0, sigma = 0.5, lambda = 0.5, intervals = 1:100, seed=123)

sc_one_series_obs <- c(5, 25, 50, 100)

sc_one_series_pars <- c(lambda=.2, sigma=.2, mu=0)

# Samples
sc_one_series_samples <- lapply(sc_one_series_obs, function(x) {
  sampling(hierarchical_student_t_oup, generate_student_set(1, student_df = 7, mu = sc_one_series_pars["mu"], sigma = sc_one_series_pars["sigma"], lambda = sc_one_series_pars["lambda"], intervals = 1:x, seed=531), chains=2, iter=2000, init=1)
}) %>% set_names(sc_one_series_obs %>% as.character)

# save(sc_one_series_samples, file="sc_one_series_samples")
# load("sc_one_series_samples")



# ('#fdd49e','#fc8d59','#d7301f','#7f0000', 'black')

 ## Results
sc_one_series_pos_seq_plot_lambda <- plot_pos_sequence(sc_one_series_samples, par="lambda") + guides(color=guide_legend(title="Sample size")) + geom_vline(xintercept=0.5, linetype="dashed") + scale_x_continuous(limits=c(0, 1)) + xlab(TeX("$\\lambda")) + scale_color_manual(values=c('#fc8d59','#d7301f','#7f0000', 'black'))

sc_one_series_pos_seq_plot_mu <- plot_pos_sequence(sc_one_series_samples, par="mu") + guides(color=guide_legend(title="Sample size")) + geom_vline(xintercept=0, linetype="dashed") + scale_x_continuous(limits=c(-.5, .5)) + xlab(TeX("$\\mu")) + scale_color_manual(values=c('#fc8d59','#d7301f','#7f0000', 'black'))

sc_one_series_pos_seq_plot_sigma <- plot_pos_sequence(sc_one_series_samples, par="sigma") + guides(color=guide_legend(title="Sample size")) + geom_vline(xintercept=0.2, linetype="dashed") + scale_x_continuous(limits=c(0, 1)) + xlab(TeX("$\\sigma")) + scale_color_manual(values=c('#fc8d59','#d7301f','#7f0000', 'black'))

sc_one_series_pos_plot_lambda <-  plot_pos(sc_one_series_samples[[sc_one_series_obs %>% max %>% as.character()]], par="lambda") + geom_vline(xintercept=0.5, linetype="dashed") + xlab(TeX("$\\lambda")) + scale_x_continuous(limits=c(0, 1))

sc_one_series_pos_plot_mu <-  plot_pos(sc_one_series_samples[[sc_one_series_obs %>% max %>% as.character()]], par="mu") + geom_vline(xintercept=0, linetype="dashed") + xlab(TeX("$\\mu"))

sc_one_series_pos_plot_sigma <-  plot_pos(sc_one_series_samples[[sc_one_series_obs %>% max %>% as.character()]], par = "sigma") + geom_vline(xintercept=0.2, linetype="dashed") + xlab(TeX("$\\sigma"))+ scale_x_continuous(limits=c(0, 1))




#### Hierarchical ####

# sc_many_n_series <- c(1, 5, 25, 50, 100, 150)
sc_many_n_series <- c(1, 5, 25, 50, 100, 150)
sc_many_series_samples <- list()
for(i in sc_many_n_series) {
  
  sc_many_series_samples[[i %>% as.character()]] <- sampling(hierarchical_student_t_oup, generate_student_set(i, student_df = 7, mu = rep(0, i), sigma = rep(0.5, i), lambda = rep(0.5,i), intervals = 1:5, seed=123), init=1, chains=2, iter=2000)
  
}

# save(sc_many_series_samples, file="sc_many_series_samples")

sc_many_series_samples <- sc_many_series_samples %>% set_names(sc_many_n_series %>% as.character)


sc_many_series_pos_plot <- plot_pos_sequence(sc_many_series_samples, par = "mu[1]") + guides(color=guide_legend(title="Number of series")) + geom_vline(xintercept = 0.5, linetype="dashed") + scale_color_brewer(type="seq" , palette = 2)
