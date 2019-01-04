library(rstan)
library(mvtnorm)
library(tidyverse)
library(invgamma)
library(reshape2)
library(magrittr)

source("OU.functions.R")


chains <- 2
iter <- 2000


options(mc.cores = parallel::detectCores())