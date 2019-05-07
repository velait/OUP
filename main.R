library(rstan)
library(mvtnorm)
library(tidyverse)
library(invgamma)
library(reshape2)
library(magrittr)
library(cowplot)
library(ggsci)
library(shinystan)
library(KScorrect)
library(Matrix)
library(matrixcalc)
library(GGally)
library(Sim.DiffProc)
library(cusp)

source("OU.functions.R")


chains <- 1
iter <- 2000


options(mc.cores = parallel::detectCores())