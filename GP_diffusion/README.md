## Literature

## See intro section in https://arxiv.org/pdf/1704.04375.pdf to get the idea behid the project.  

## This should have an approachable overview of stochastic differential equations: https://link.springer.com/chapter/10.1007%2F978-3-030-01768-2_16

## Gaussian process 
# Good introductory material:
# - https://distill.pub/2019/visual-exploration-gaussian-processes/
# - http://www.tmpl.fi/gp/
# - GPs in Stan
#    - https://mc-stan.org/docs/2_19/stan-users-guide/gaussian-process-regression.html
#    - https://betanalpha.github.io/assets/case_studies/gp_part1/part1.html
# - https://en.wikipedia.org/wiki/Gaussian_process
# - THE GP book: http://gaussianprocess.org/gpml/chapters/RW.pdf



## Scripts 

# setup.R
Get the necessary R packages

# hitchip_functions.R
Functions. Probably contains some redundant material.

# GP_example.R
An independent and minimal example of non-parametric regression with Gaussian processes.

# GP.stan
Stan program for fitting a Gaussian process

# HitChip_data.R
Wrangles data in HitChip atlas phyloseq, output is abundance and meta data tables, written in .txt files. 