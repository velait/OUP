# Fix "ground state:"  at grazing rate c = 1 --> oup parameters --> test predictability


set.seed(1)
times <- seq(from = 0, to = 1000, by = 1)
N_series <- 100


# Ground state data

r <- 1
K <- 10
cs <- 1
h <- 1
sigma <- 0.07


y <- rep(NA, length(times))
y[1] <- 8

cs <- 1 + (2.6771 - 1)*(times/max(times))


ground_state <- ews_generator(8, times, c = 1, milstein = T)
  






## Fit OUP to ground state

# frequentist
ou.lik <- function(x) {
  function(theta1,theta2,theta3) {
    n <- length(x)
    dt <- deltat(x)
    # dt <- .1
    -sum(dcOU(x=x[2:n], Dt=dt, x0=x[1:(n-1)],
              theta=c(theta1,theta2,theta3), log=TRUE))
  }
}
ground_oup_pars <- lapply(1, function(j) {
  
  ou.fit <- mle(ou.lik(ground_state),
                start=list(theta1=1,theta2=0.5,theta3=0.2),
                method="L-BFGS-B",lower=c(0,1e-5,1e-3), upper=c(1,1,1))
  ou.coe <- coef(ou.fit)
  coefs <- c(ou.coe[2], ou.coe[1]/ou.coe[2], ou.coe[3]) %>% set_names(c("lambda", "mu", "sigma"))
  
  
  # # Next time point likelihood
  # next_m <- coefs["mu"] - (coefs["mu"] - my_data[j + window_length - 1])*exp(-coefs["lambda"])
  # next_sd <- ((coefs["sigma"])/(2*coefs["sigma"]))*(1-exp(-2*coefs["lambda"]))
  # LL <- dnorm(x,
  #             mean = next_m,
  #             sd = next_sd) 
  # 
  # 
  
})[[1]]

# Stan
oup_model <- stan_model("stan_models/oup_single_transition.stan")
oup_samples <- sampling(oup_model, list(T = length(times),
                                        time = times,
                                        Y = ground_state),
                        chains = 2, iter = 1000)

oup_pars <- c("lambda", "mu", "sigma")
oup_res <- lapply(oup_pars, function(p) {
  
  get_stan_results(oup_samples, p)
  
}) %>% do.call(rbind, .) %>% pull(mode) %>% set_names(oup_pars)



## PPC: Compare data to inferred parameters
data.frame(ground_state[1:1000], oup_generator(N = 1000, times = 1:1000, y0 = 9, mu = oup_res["mu"], lambda = oup_res["lambda"], sigma = oup_res["sigma"])) %>%
  set_colnames(c("ground", "infer")) %>% 
  mutate(x = 1:1000) %>% 
  melt(id.vars = "x") %>% 
  ggplot(aes(x = x, y =  value, color = variable)) + 
    geom_line()+ coord_cartesian(xlim = 1:100)



ind <- 1
series <- ews_set[, ind]

# coefs <- ground_oup_pars
coefs <- c(1.46, 8.76, 1.42) %>% set_names(c("lambda", "mu", "sigma"))

ground_state_comparison <- lapply(2:length(series), function(i) {
  
  prev <- series[i-1]
  present <- series[i]
  
  next_m <- coefs["mu"] - (coefs["mu"] - prev)*exp(-coefs["lambda"])
  next_sd <- ((coefs["sigma"]^2)/(2*coefs["lambda"]))*(1-exp(-2*coefs["lambda"]))
  
  # Likelihood
  LL <- dnorm(present,
              mean = next_m,
              sd = next_sd)
  
  # Which interquartile range is the observation?
  Q <- abs(0.5 - pnorm(present, mean = next_m, sd = next_sd))*2
  
  
  c(likelihood = LL, iqr = Q, timeindex = times[i], upper = next_m + next_sd, lower = next_m - next_sd , true = present)
  
}) %>% do.call(rbind, .)


# ground_state_comparison[, "likelihood"] <- (ground_state_comparison[, "likelihood"] - mean(ground_state_comparison[, "likelihood"]))/sd(ground_state_comparison[, "likelihood"])
# 
# ground_state_comparison[, "likelihood"] <- -ground_state_comparison[, "likelihood"]


ground_state_comparison %>% 
  as.data.frame() %>% 
  # filter(timeindex < 900) %>%
  ggplot() + 
  geom_line(aes(x = timeindex, y = true)) + 
  geom_errorbar(aes(x = timeindex, ymin = lower.mu, ymax = upper.mu), color = "royalblue", alpha = 0.25) +
  geom_line(aes(x = timeindex, y = likelihood), color = "red") +
  geom_line(aes(x = timeindex, y = iqr), color = "green") +
  labs(subtitle = "Red = likelihood of an observation given the previous one; \n Blue = mean +- 1 sd of the transition density: \n Green = Smallest IQR containing the next observation")
  # coord_cartesian(ylim = 0:10)
  






