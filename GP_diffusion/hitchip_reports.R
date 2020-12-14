# Test GP interpolation approach on HITCHIP Atlas

data("atlas1006")
hitchip_pseq <- get_edit_hitchip_pseq(atlas1006)
hitchip_data <- get_hitchip_data(hitchip_pseq, pcoa = FALSE, transformation = "log10")

# hitchip_otu <- "Akkermansia"

n_samples




core_taxa <- core(hitchip_pseq %>% transform("compositional"), detection = .05, prevalence = .1) %>% 
  taxa_names
bimodal_taxa <- c("Anaerofustis" = 0.4, 
                   "Aneurinibacillus" = 0.3, 
                   "Aquabacterium" = 0.5, 
                   "Bacteroides_intestinalis_et_rel." = .5, 
                   "Bacteroides_fragilis_et_rel." = 3.1,
                   "Burkholderia" = 0.5, 
                   "Dialister" = 3.25, 
                   "Leminorella" = 0.5, 
                   "Prevotella_melaninogenica_et_rel." = 4.5, 
                   "Prevotella_oralis_et_rel." = 4, 
                   "Serratia" = 0.5, 
                   "Uncultured_Bacteroidetes" = 0.5, 
                   "Uncultured_Clostridiales_II" = 3.4, 
                   "Uncultured_Selenomonadaceae" = 0.5, 
                   "Wissella_et_rel." = 0.45) %>% names


test_taxa <- c(core_taxa, bimodal_taxa)
test_taxa <- test_taxa[!duplicated(test_taxa)]

stationary_density_plots <- lapply(test_taxa, function(taxon) {
  print(taxon)
  
  p <- hitchip_stationary_density_plotter(hitchip_data, taxon)
  
  return(p)
}) %>%
  set_names(test_taxa)

