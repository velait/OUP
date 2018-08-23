# Example for the IDA paper

#### SETTIGS ####
library(rstan)
library(tidyverse)
library(magrittr)
library(MCMCpack)
library(mvtnorm)
library(reshape2)
library(latex2exp)
library(cowplot)


theme_set(theme_bw(12))
options(mc.cores = parallel::detectCores())
set.seed(11235)

source("OU.functions.R")


iter <- 2000
chains <- 2

#### EXAMPLE ####
## Data
example_data <- generate_student_set(1, student_df = 7, mu = 0, sigma = 0.5, lambda = 0.5, intervals = 1:150, seed=69)


# Model
pooled_student_t_oup <- stan_model("pooled_student_t_oup.stan")

# Plot
example_plot <- ggplot() + geom_line(data=example_data[["Y"]] %>% t %>%  as_tibble, aes(y=V1, x=1:150)) + labs(x="Time", y="") + scale_y_continuous(limits = c(-3, 3))

# Sample
example_samples <- sampling(pooled_student_t_oup, example_data, chains=2, iter=2000)

#### Results

# Samples
example_pos <- sapply(c("lambda", "mu", "sigma"), function(x) rstan::extract(example_samples, x)) %>% set_names(c("lambda", "mu", "sigma"))

# melt_pos <- do.call(cbind, example_pos) %>% as_tibble() %>% set_colnames(c("lambda", "mu", "sigma")) %>% melt

# collect to a df
df <- do.call(cbind, example_pos) %>% as_tibble() %>% set_colnames(c("lambda", "mu", "sigma"))

# list for plots
example_pos_plot <- list()

example_pos_plot[["mu"]] <- ggplot(df, aes(x=mu)) + stat_density(geom="line") + labs(y="", x=expression(~mu)) + scale_x_continuous(limits=c(-.5, .5)) + geom_vline(xintercept = 0, linetype="dashed")

example_pos_plot[["lambda"]] <- ggplot(df, aes(x=lambda)) + stat_density(geom="line") + labs(y="Density", x=expression(~lambda)) + scale_x_continuous(limits=c(0, 1.1), breaks=c(0, .25, .5, .75, 1)) + geom_vline(xintercept = 0.5, linetype="dashed")

example_pos_plot[["sigma"]] <- ggplot(df, aes(x=sigma)) + stat_density(geom="line") + labs(y="", x=expression(~sigma)) + scale_x_continuous(limits=c(0, 1.2)) + geom_vline(xintercept = 0.5, linetype="dashed")


# Arrange to a grid
posterior_row <- plot_grid(example_pos_plot[["lambda"]], example_pos_plot[["mu"]], example_pos_plot[["sigma"]], nrow=1)


example_plot <- plot_grid(example_plot, posterior_row, labels = c("A", "B"), nrow = 2, rel_heights = c(1, 1))


# example_plot %>% print



