R_packages <- c("magrittr", "tidyverse", "reshape2", "rstan")

# Install missing packages
missing_packages <- R_packages[!(R_packages %in% installed.packages()[,"Package"])]
if(length(missing_packages) != 0) install.packages(missing_packages)


