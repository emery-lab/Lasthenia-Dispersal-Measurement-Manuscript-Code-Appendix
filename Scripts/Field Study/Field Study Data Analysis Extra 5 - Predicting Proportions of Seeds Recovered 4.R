#### Field Dispersal Data Analysis Extra 5: Predicting Proportions of Seed Recovered 4

# Courtney Van Den Elzen

# Initiated: October 2019
# Latest Update: January 2023 (syntax update for Jan 2023 package version compatibility)

# Description: Summarize the total numbers of seeds found in the field for each 
#              maternal plant and calculate predicted proportions


# WARNING: You must first run the Field Study Data Analysis Extra 3 & 4 scripts

### ------------------
### Read In Maternal Dispersal Kernel Data -----

kernel_seeds <- read_csv("./Data/Field Study/Field Study Raw Data - Dispersal Kernels.csv",
                         col_types = cols(Transect = col_character(), 
                                          Mom = col_character())) %>% 
  dplyr::select(-c(Distance, Angle_Comp, Angle_Deg_fromN, Angle_Rad_fromN, 
                   Angle_Deg_fromE, Angle_Rad_fromE, Notes)) 


### ------------------
### Clean Data -----

# Drop any plants that I did not do collections on
kernel_seeds_good <- kernel_seeds %>% tidyr::drop_na(Pin)

# Drop any plants with fewer than 3 pins measured
kernel_seeds_total <- 
     dplyr::summarize(group_by(kernel_seeds_good, `Full ID`), 
                      n_pins = n()) %>% 
     dplyr::filter(n_pins > 3)

### ------------------
### Calculate Predicted Proportions of Recovered Seeds -----

# Calculate predicted proportions of viable seeds recovered
pred_w_recovered <- 
     dplyr::right_join(mat_inflor_data, 
                       kernel_seeds_total, 
                       by = c("Mat_ID" = "Full ID")) %>% 
     mutate(pois_prop_pred = n_pins/pois_nseed_pred, 
            pois_prop_pred_lower = n_pins/pois_nseed_pred_lower, 
            pois_prop_pred_upper = n_pins/pois_nseed_pred_upper)

# Break down by species
pred_w_recovered_cal <- pred_w_recovered %>% filter(Species == "cal")
pred_w_recovered_fre <- pred_w_recovered %>% filter(Species == "fre")
pred_w_recovered_gla <- pred_w_recovered %>% filter(Species == "gla")

### ------------------
### Plot Predicted Proportions -----

# Plot the predicted proportions of viable seeds recovered with prediction interval bars around.
ggplot(pred_w_recovered, aes(x=Diameter, y=pois_prop_pred, colour = Species)) +
  geom_pointrange(aes(ymin = pois_prop_pred_lower, 
                      ymax=pois_prop_pred_upper), 
                  position=position_jitter(width=0.5)) 

### ------------------
### Predicted Proportions (Means and Standard Deviations) by Species -----
# Means and standard deviations of predictions by species
 
# L. californica
mean(pred_w_recovered_cal$pois_prop_pred)
sd(pred_w_recovered_cal$pois_prop_pred)

# L. fremontii
mean(pred_w_recovered_fre$pois_prop_pred)
sd(pred_w_recovered_fre$pois_prop_pred)

# L. glaberrima
mean(pred_w_recovered_gla$pois_prop_pred)
sd(pred_w_recovered_gla$pois_prop_pred)

