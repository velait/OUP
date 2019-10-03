# Test non-stationary cusp catastrophy model

grid <- seq(from = 0, to = 100, by = 0.01)

r <- rep(1, length(grid))
alpha <- seq(from = -3, to = 3, length.out = length(grid))
beta <- rep(1.5, length(grid))
lambda <- rep(0, length(grid))
epsilon <- rep(.01, length(grid))


inits <- c(0)

df <- lapply(inits, function(y0) non_stat_shoji_generator(y0 = y0, times = grid, r, alpha, beta, lambda, epsilon, seed = NULL)) %>%
  do.call(cbind, .) %>%
  as.data.frame() %>%
  cbind(time = grid) %>%
  set_colnames(c("y", "x"))
# set_colnames(c("y1", "y2", "x"))

# thin_df <- thin(df, modulo = .1)

thin_df <- df

diff_plot <- thin_df %>%
  melt(id.vars = "x") %>%
  ggplot(aes(x = x, y = value, color = variable)) +
  geom_line() +
  # scale_color_tron() +
  labs(title = "", x = "Time", y = "Abundance") +
  theme(legend.position = "none")

diff_plot