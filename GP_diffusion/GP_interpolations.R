# Test GP interpolation approach on different data sets.

## 1D HitChip **************************** ####

data("atlas1006")
hitchip_pseq <- get_edit_hitchip_pseq(atlas1006)
hitchip_data <- get_hitchip_data(hitchip_pseq, pcoa = FALSE)

hitchip_otu <- "Dialister"

length_scale_par <- 1
variance_par <- .25
epsilon <- .25


x_pred <- range_seq(-5:0, length_out = 50)


mean_c <- mean(hitchip_data$abundances[, hitchip_otu])

## One dimensional
one_dimensional_hitchip_drift <- GP_imputation(X = hitchip_data$delta_df %>% 
                                                 filter(otu == hitchip_otu) %>% 
                                                 pull(x), 
                                               Y = hitchip_data$delta_df %>% 
                                                 filter(otu == hitchip_otu) %>% 
                                                 pull(unit_dx),
                                               X_pred = x_pred,
                                               rho = length_scale_par,
                                               sigma = variance_par,
                                               epsilon, 
                                               mean_c = mean_c,
                                               add_linear = TRUE)

one_dimensional_hitchip_diffusion <- GP_imputation(X = hitchip_data$delta_df %>% 
                                                 filter(otu == hitchip_otu) %>% 
                                                 pull(x), 
                                               Y = hitchip_data$delta_df %>% 
                                                 filter(otu == hitchip_otu) %>% 
                                                 add_difference_sq_1D() %>% 
                                                 pull(difference_sq)
                                                 ,
                                               X_pred = x_pred,
                                               rho = length_scale_par,
                                               sigma = variance_par,
                                               epsilon, 
                                               mean_c, add_linear = TRUE)


hitchip_1d_plots <- GP_1D_diffusion_plotter(GP_drift = one_dimensional_hitchip_drift, 
                                            GP_diffusion = one_dimensional_hitchip_diffusion, 
                                            full_data = hitchip_data$abundances[, hitchip_otu], 
                                            variable = hitchip_otu)





## 2D HitChip **************************** ####

hitchip_otu1 <- "Bacteroides_fragilis_et_rel."
hitchip_otu2 <- "Akkermansia"

# hitchip_otu1 <- "Axis.1"
# hitchip_otu2 <- "Axis.2"

hitchip_pairwise_df <- get_pairwise_df(hitchip_data$delta_df,
                                       x_otu = hitchip_otu1, 
                                       y_otu = hitchip_otu2)


# normalize_pairwise_df <- function(df) {
#   
#   
#   x_location <- min(df$x)
#   x_scale <- max((df$x - min(df$x)))
#   
#   df$unit_dx <- (df$unit_dx - x_location)/x_scale
#   df$x <- (df$x - x_location)/x_scale
#   
#   
#   
#   
#   y_location <- min(df$y)
#   y_scale <- max((df$y - min(df$y)))
#   
#   df$unit_dy <- (df$unit_dy - y_location)/y_scale
#   df$y <- (df$y - y_location)/y_scale
#   
#   
#   return(df)
# }
# hitchip_pairwise_df <- normalize_pairwise_df(hitchip_pairwise_df)

resolution <- 10

full_data <- hitchip_data$abundances[, c(hitchip_otu1, hitchip_otu2)]


# Prediction points
hitchip_pred_range <- expand.grid(x = range_seq(range(full_data[, 1]) + c(-.25, 25), length_out = resolution),
                                  y = range_seq(range(full_data[, 2])+ c(-.25, 25), length_out = resolution))


mean_c <- hitchip_pred_range %>% colMeans()

## Multi-dimensional
mv_hitchip_drift <- GP_imputation(X = hitchip_pairwise_df[, c("x", "y")], 
                                  Y = hitchip_pairwise_df[, c("unit_dx", "unit_dy")],
                                  X_pred = hitchip_pred_range,
                                  rho = c(.5, .5), sigma = .1, epsilon, 
                                  add_linear = TRUE,  mean_c = mean_c)



# p <- full_join(grid_average(hitchip_pairwise_df[, c("x")],
#                  hitchip_pairwise_df[, c("y")],
#                  hitchip_pairwise_df[, c("unit_dx")]),
# grid_average(hitchip_pairwise_df[, c("x")],
#              hitchip_pairwise_df[, c("y")],
#              hitchip_pairwise_df[, c("unit_dy")]), by = c("x", "y")) %>%
#   ggplot()  +
#   geom_segment(aes(x = x, y = y, xend = x + z.x, yend = y + z.y),
#                arrow = arrow(length = unit(0.2, "cm"), ends="last"), size = .5)



mv_hitchip_diffusion <- GP_imputation(X = hitchip_pairwise_df[, c("x", "y")], 
                                      Y = hitchip_pairwise_df %>% 
                                        add_difference_sq_2D() %>% 
                                        pull(difference_sq),
                                      X_pred = hitchip_pred_range,
                                      rho = .5,
                                      sigma = .1,
                                      epsilon)




hitchip_mv_plots <- GP_2D_diffusion_plotter(GP_drift = mv_hitchip_drift, 
                                            GP_diffusion = mv_hitchip_diffusion, 
                                            full_data = full_data, 
                                            variables = c(hitchip_otu1, 
                                                          hitchip_otu2),
                                            arrow_scaling = 1)



plot_grid(hitchip_mv_plots[["drift_diffusion"]] + labs(subtitle = paste0("HitChip; x = ", hitchip_otu1, "; \n y = ", hitchip_otu2)),
          hitchip_mv_plots[["drift_variance"]]+ labs(subtitle = "HitChip; uncertainty in drift"),
          hitchip_mv_plots[["diffusion_variance"]]+ labs(subtitle = "HitChip; uncertainty in diffusion"), ncol = 3)


plot_grid(hitchip_mv_plots[["drift_diffusion"]] + coord_cartesian(xlim = range(full_data[, 1]), 
                                                                  ylim = range(full_data[, 2])),
          hitchip_mv_plots[["base"]] + coord_cartesian(xlim = range(full_data[, 1]), 
                                                       ylim = range(full_data[, 2])))


## 1D David ****************************** ####

david_data <- get_david_data("log10")
david_otu <- "Dialister"

## One dimensional
one_dimensional_david_drift <- GP_imputation(X = rbind(david_data$delta_df$A, 
                                                       david_data$delta_df$B) %>% 
                                                 filter(otu == david_otu) %>% 
                                                 pull(x), 
                                               Y = rbind(david_data$delta_df$A, 
                                                         david_data$delta_df$B) %>% 
                                                 filter(otu == david_otu) %>% 
                                                 pull(unit_dx),
                                               X_pred = range_seq(rbind(david_data$delta_df$A, 
                                                                        david_data$delta_df$B) %>% 
                                                                    filter(otu == david_otu) %>% 
                                                                    pull(x), length_out = 50),
                                               rho = length_scale_par, sigma = variance_par, epsilon)

one_dimensional_david_diffusion <- GP_imputation(X = rbind(david_data$delta_df$A, 
                                                           david_data$delta_df$B) %>% 
                                                     filter(otu == david_otu) %>% 
                                                     pull(x), 
                                                   Y = rbind(david_data$delta_df$A, 
                                                             david_data$delta_df$B) %>% 
                                                     filter(otu == david_otu) %>% 
                                                     add_difference_sq_1D() %>% 
                                                     pull(difference_sq),
                                                   X_pred = range_seq(rbind(david_data$delta_df$A, 
                                                                            david_data$delta_df$B) %>% 
                                                                        filter(otu == david_otu) %>% 
                                                                        pull(x), length_out = 50),
                                                   rho = length_scale_par, sigma = variance_par, epsilon)


david_1d_plots <- GP_1D_diffusion_plotter(GP_drift = one_dimensional_david_drift, 
                                            GP_diffusion = one_dimensional_david_diffusion, 
                                            full_data = cbind(one_dimensional_david_drift$training$X, 
                                                              one_dimensional_david_drift$training$Y), 
                                            variable = david_otu)



## 2D David ****************************** ####

david_otu1 <- "Bacteroides"
david_otu2 <- "Prevotella"

david_pairwise_df <- get_pairwise_df(david_data$delta_df$A,
                                       x_otu = david_otu1, 
                                       y_otu = david_otu2)

resolution <- 20

# Prediction points
david_pred_range <- expand.grid(x = range_seq(david_pairwise_df[, "x"], length_out = resolution), 
                                  y = range_seq(david_pairwise_df[, "y"], length_out = resolution))

# hitchip_pred_range <- expand.grid(x = seq(from = -0.3, to = .5, length.out = resolution), 
#                                   x = seq(from = -0.4, to = .5, length.out = resolution))

full_data <- david_data$adundances[, c(david_otu1, david_otu2)]

## Multi-dimensional
mv_david_drift <- GP_imputation(X = david_pairwise_df[, c("x", "y")], 
                                  Y = david_pairwise_df[, c("unit_dx", "unit_dy")],
                                  X_pred = david_pred_range,
                                  rho = c(1, 1), sigma = .25, epsilon)

mv_david_diffusion <- GP_imputation(X = david_pairwise_df[, c("x", "y")], 
                                      Y = david_pairwise_df %>% 
                                        add_difference_sq_2D() %>% 
                                        pull(difference_sq),
                                      X_pred = david_pred_range,
                                      rho = length_scale_par, sigma = variance_par, epsilon)


david_mv_plots <- GP_2D_diffusion_plotter(GP_drift = mv_david_drift, 
                                            GP_diffusion = mv_david_diffusion, 
                                            full_data = full_data, 
                                            variables = c(david_otu1, 
                                                          david_otu2),
                                            arrow_scaling = 1)


# david_mv_plots$drift_diffusion
# david_mv_plots$diffusion_variance


## 1D Ravel ****************************** ####

ravel_data <- get_ravel_data()

ravel_otu <- "Bifidobacterium"
ravel_subject <- 1:32

## One dimensional
one_dimensional_ravel_drift <- GP_imputation(X = ravel_data$delta_df %>% 
                                               filter(otu == ravel_otu, subject %in% ravel_subject) %>% 
                                               pull(x), 
                                             Y =  ravel_data$delta_df %>% 
                                               filter(otu == ravel_otu, subject %in% ravel_subject) %>% 
                                               pull(unit_dx),
                                             X_pred = range_seq( ravel_data$delta_df %>% 
                                                                  filter(otu == ravel_otu, subject %in% ravel_subject) %>% 
                                                                  pull(x), length_out = 50),
                                             rho = length_scale_par, sigma = variance_par, epsilon)

one_dimensional_ravel_diffusion <- GP_imputation(X = ravel_data$delta_df %>% 
                                                   filter(otu == ravel_otu, subject %in% ravel_subject) %>% 
                                                   pull(x), 
                                                 Y = ravel_data$delta_df %>% 
                                                   filter(otu == ravel_otu, subject %in% ravel_subject) %>% 
                                                   add_difference_sq_1D() %>% 
                                                   pull(difference_sq),
                                                 X_pred = range_seq(ravel_data$delta_df %>% 
                                                                      filter(otu == ravel_otu, subject %in% ravel_subject) %>% 
                                                                      pull(x), length_out = 50),
                                                 rho = length_scale_par, sigma = variance_par, epsilon)


ravel_1d_plots <- GP_1D_diffusion_plotter(GP_drift = one_dimensional_ravel_drift, 
                                          GP_diffusion = one_dimensional_ravel_diffusion, 
                                          full_data = ravel_data$abundances[, ravel_otu], 
                                          variable = ravel_otu, 
                                          scaling = 30)



## Loop over all ravel subjects

ravel_loop <- lapply(1:32, function(i) {
  
  ravel_subject <- i
  
  ## One dimensional
  one_dimensional_ravel_drift <- GP_imputation(X = ravel_data$delta_df %>% 
                                                 filter(otu == ravel_otu, subject %in% ravel_subject) %>% 
                                                 pull(x), 
                                               Y =  ravel_data$delta_df %>% 
                                                 filter(otu == ravel_otu, subject %in% ravel_subject) %>% 
                                                 pull(unit_dx),
                                               X_pred = range_seq( ravel_data$delta_df %>% 
                                                                     filter(otu == ravel_otu, subject %in% ravel_subject) %>% 
                                                                     pull(x), length_out = 50),
                                               rho = length_scale_par, sigma = variance_par, epsilon)
  
  one_dimensional_ravel_diffusion <- GP_imputation(X = ravel_data$delta_df %>% 
                                                     filter(otu == ravel_otu, subject %in% ravel_subject) %>% 
                                                     pull(x), 
                                                   Y = ravel_data$delta_df %>% 
                                                     filter(otu == ravel_otu, subject %in% ravel_subject) %>% 
                                                     add_difference_sq_1D() %>% 
                                                     pull(difference_sq),
                                                   X_pred = range_seq(ravel_data$delta_df %>% 
                                                                        filter(otu == ravel_otu, subject %in% ravel_subject) %>% 
                                                                        pull(x), length_out = 50),
                                                   rho = length_scale_par, sigma = variance_par, epsilon)
  
  
  return(list(drift = one_dimensional_ravel_drift, 
              diffusion = one_dimensional_ravel_diffusion))
  
})


loop_plot <- GP_1D_diffusion_plotter_set(ravel_loop,
                            full_data = ravel_data$abundances[, ravel_otu], 
                            scaling = 10, 
                            variable = ravel_otu)


## 2D Ravel ****************************** ####


ravel_otu1 <- "Bacteroides"
ravel_otu2 <- "Prevotella"

ravel_pairwise_df <- get_pairwise_df(ravel_data$delta_df,
                                     x_otu = ravel_otu1, 
                                     y_otu = ravel_otu2)

resolution <- 20

# Prediction points
ravel_pred_range <- expand.grid(x = range_seq(ravel_pairwise_df[, "x"], length_out = resolution), 
                                y = range_seq(ravel_pairwise_df[, "y"], length_out = resolution))

# hitchip_pred_range <- expand.grid(x = seq(from = -0.3, to = .5, length.out = resolution), 
#                                   x = seq(from = -0.4, to = .5, length.out = resolution))

full_data <- ravel_data$abundances[, c(ravel_otu1, ravel_otu2)]

ravel_range <- 1:100

## Multi-dimensional
mv_ravel_drift <- GP_imputation(X = ravel_pairwise_df[ravel_range, c("x", "y")], 
                                Y = ravel_pairwise_df[ravel_range, c("unit_dx", "unit_dy")],
                                X_pred = ravel_pred_range,
                                rho = c(1, 1), sigma = .25, epsilon)

mv_ravel_diffusion <- GP_imputation(X = ravel_pairwise_df[ravel_range, c("x", "y")], 
                                    Y = ravel_pairwise_df[ravel_range, ] %>% 
                                      add_difference_sq_2D() %>% 
                                      pull(difference_sq),
                                    X_pred = ravel_pred_range,
                                    rho = length_scale_par, sigma = variance_par, epsilon)


ravel_mv_plots <- GP_2D_diffusion_plotter(GP_drift = mv_ravel_drift, 
                                          GP_diffusion = mv_ravel_diffusion, 
                                          full_data = full_data, 
                                          variables = c(ravel_otu1, 
                                                        ravel_otu2),
                                          arrow_scaling = 1)


# david_mv_plots$drift_diffusion
# david_mv_plots$diffusion_variance
