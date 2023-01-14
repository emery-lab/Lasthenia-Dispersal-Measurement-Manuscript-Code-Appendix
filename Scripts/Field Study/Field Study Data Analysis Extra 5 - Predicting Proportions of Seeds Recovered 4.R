#' # Field Study 2019 Jepson Prairie
#' # Recovered Seed Proportions Part 3
#' 
#' ### Part 4/4: Summarize the total numbers of seeds found in the field for each maternal plant and calculate predicted proportions
library(tidyverse)
#' Courtney Van Den Elzen
#' 
#' October 21st 2019 - Nov 22nd 2019
#' 
#' Part 3 must be run first

#' #### Read in the maternal dispersal kernel data
kernel_seeds <- read_csv("/Users/Courtney/Documents/Thesis Chapter 2/Field Study/2019 Full Study Files/Data/2019 Field Data Jepson Dispersal Kernels.csv",
                         col_types = cols(Transect = col_character(), Mom = col_character())) %>% 
  dplyr::select(-c(Distance, Angle_Comp, Angle_Deg_fromN, Angle_Rad_fromN, Angle_Deg_fromE, Angle_Rad_fromE, Notes)) 

#' Drop any plants that I did not do collections on
kernel_seeds_good <- kernel_seeds %>% tidyr::drop_na(Pin)

#' Drop any plants with fewer than 3 pins measured
kernel_seeds_total <- dplyr::summarize(group_by(kernel_seeds_good, `Full ID`), n_pins = n()) %>% filter(n_pins > 3)

#' Calculate predicted proportionsof viable seeds recovered
pred_w_recovered <- right_join(mat_inflor_data, kernel_seeds_total, by = c("Mat_ID" = "Full ID")) %>% mutate(pois_prop_pred = n_pins/pois_nseed_pred, 
                                                                                                             pois_prop_pred_lower = n_pins/pois_nseed_pred_lower, 
                                                                                                             pois_prop_pred_upper = n_pins/pois_nseed_pred_upper)
pred_w_recovered_cal <- pred_w_recovered %>% filter(Species == "cal")
pred_w_recovered_fre <- pred_w_recovered %>% filter(Species == "fre")
pred_w_recovered_gla <- pred_w_recovered %>% filter(Species == "gla")

#' Plot the predicted proportions of viable seeds recovered with prediction interval bars around.
ggplot(pred_w_recovered, aes(x=Size, y=pois_prop_pred, colour = Species)) +
  geom_pointrange(aes(ymin =pois_prop_pred_lower, ymax=pois_prop_pred_upper), position=position_jitter(width=0.5)) 

#' Means and standard deviations of predictions by species
#' 
#' *L. californica*
mean(pred_w_recovered_cal$pois_prop_pred)
sd(pred_w_recovered_cal$pois_prop_pred)

#' *L. fremontii*
mean(pred_w_recovered_fre$pois_prop_pred)
sd(pred_w_recovered_fre$pois_prop_pred)

#' *L. glaberrima*
mean(pred_w_recovered_gla$pois_prop_pred)
sd(pred_w_recovered_gla$pois_prop_pred)

#' For knitting to html: knitr::spin("/Users/Courtney/Documents/Thesis Chapter 1/Field Study/2019 Full Study Files/Scripts/Predicting Found Proportions 4.R")

