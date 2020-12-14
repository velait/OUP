
# ********************************* #

# In https://arxiv.org/pdf/1704.04375.pdf it is shown (old result, however)
# that the drift and diffusion coefficients can be approximated with the difference between consecutive observations.
# Let (x_i, y_i) be observations of a dynamical system. Then unit_dx = (x_{i+1} - x_{i})/(y_{i+1} - y_{i})
# approximates the drift at x_{i} and unit_dx^2 approximates the diffusion.

# In HitChip Atlas we have some longitudinal observations from some individuals.
# In order to learn the SDE that governs the dynamics of a bacterial species on population level
# we can combine all the available abundance changes and then do non-parametric regression on the differences.
# So far I have only done it separately for the drift and diffusion coefficient but I guess the learning could be done simulatneously too. 

# Functions in this script produce data from HitChip Atlas phyloseq objects that allows the regression described above.
# See hitchip_functions.R for details. 

# ********************************* #


data("atlas1006")

# Get phyloseq and remove duplicated observations
hitchip_pseq <- get_edit_hitchip_pseq(atlas1006)

# Get subjects with multiple observations
# Compute abundances differences
hitchip_data <- get_hitchip_data(hitchip_pseq,
                                 transformation = "log10",
                                 pcoa = FALSE)



# Check data. Columns:
# dt = time change
# dx = abundance change
# x = intial point
# otu = species
# subject = id (same as in meta data)
# unit_dx = dx/dt
head(hitchip_data$delta_df)



# All sample abundances. 
head(hitchip_data$abundances[, 1:5])

# All longintudinal sample abundances 
head(hitchip_data$longitudinal_abundances[, 1:5])

# Meta data
head(meta(hitchip_pseq))



## Write to .csv
write.table(hitchip_data$delta_df, file = "GP_diffusion/delta_df.txt")
write.table(hitchip_data$abundances, file = "GP_diffusion/all_abundances.txt")
write.table(hitchip_data$longitudinal_abundances, file = "GP_diffusion/longitudinal_abundances.txt")
write.table(meta(hitchip_pseq), file = "GP_diffusion/meta_data.txt")