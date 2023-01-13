#### 2019 Field Dispersal Data Analysis 4: Combining Maternal, Surrounding Vegetation, Wind Speed, and Dispersal Data

# Courtney Van Den Elzen

# Initiated: December 2019
# Latest Update: January 2023 (syntax update for Jan 2023 package version compatibility)


### ------------------
### Load Packages
library(tidyverse)

### ------------------
### Read In All Cleaned Data

# maternal inflorescence count data
midf_n_inf <- read_csv("./Data/Field Study Cleaned Data - Maternal Inflorescences.csv")

#maternal plant traits totals
mpdf_totals_means_var <- read_csv("./Data/Field Study Cleaned Data - Maternal Plant Mass Total Mean and Var Lengths.csv")

#maternal plant traits all branch lengths
mpdf_long <- read_csv("./Data/Field Study Cleaned Data - Maternal Plant Branch Lengths Long Form.csv")

#surrounding vegetation mass and length (max, mean) data
svdf_mass_len <- read_csv("./Data/Field Study Cleaned Data - Surrounding Veg Mass Max Mean Length.csv")

#handheld windspeed 
hwdf_means <- read_csv("./Data/Field Study Cleaned Data - Wind Speed Means.csv")

# dispersal (kernel) data
kernels <- read_csv("./Data/Field Study Raw Data - Dispersal Kernels.csv") %>% dplyr::select(-c(Species))
kernels_w_disp <- read_csv("./Data/Field Study Raw Data - Dispersal Kernels w Dispersion.csv") %>% dplyr::select(-c(Species))
kernels_w_disp$...1 <- NULL

### ------------------
### Combine All Cleaned Data

# Combine the maternal inflorescence counts with the maternal plant traits data
mpdf_midf <- full_join(mpdf_totals_means_var, midf_n_inf, by = c("Full ID" = "Plant_ID"))

# Look at the correlation structure within the maternal plant traits.
cor(drop_na(dplyr::select(mpdf_midf, `Mass (mg)`, `Main Stem Length (mm)`, Tot_len, Mean_len, Var_len, n_inflor)))

# Further clean and examine the surrounding vegetation data
svdf_mass_len <- 
     svdf_mass_len %>% 
     dplyr::select(Plant_ID, Species, Hub, Transect, Mom, max_length_mm, Mass_g, mean_length_mm, veg_density_kg_m2, Notes)

# Look at the correlation structure within the surrounding vegetation data
cor(drop_na(dplyr::select(svdf_mass_len, max_length_mm, Mass_g, mean_length_mm, veg_density_kg_m2)))

# Further clean and examine the handheld windspeed data
cor(drop_na(dplyr::select(hwdf_means, mean_speed_1m, mean_speed_10cm, mean_prop_speed)))

# Combine the surrounding vegetation data with the maternal plant data
sv_mp_midf <- full_join(svdf_mass_len, mpdf_midf, by = c("Plant_ID" = "Full ID", "Species", "Hub", "Transect", "Mom"))


colnames(sv_mp_midf) <- c("Plant_ID", "Species", "Hub", "Transect", "Mom", 
                          "max_surveg_length_mm", "surveg_mass_g", 
                          "mean_surveg_length_mm", "veg_density_kg_msqr", 
                          "Notes", "mat_mass_mg", "main_stem_length_mm", 
                          "tot_mat_length_mm", "mean_inflor_height_mm", 
                          "var_inflor_height_mm", "non_na_count", "n_inflor")

# Combine the windspeed data with the other datasets
hw_sv_mp_midf <- 
     hwdf_means %>% 
     full_join(sv_mp_midf, by = c(`Full ID` = "Plant_ID", "SppID" = "Species", "Hub", "Transect", "Mom")) %>% 
     dplyr::select(-Notes) %>% 
     dplyr::filter(!is.na(SppID))

# Join the kernel data in with the other data
kern_hw_sv_mp_midf <- full_join(kernels, hw_sv_mp_midf)
kern_w_disp_hw_sv_mp_midf <- full_join(kernels_w_disp, hw_sv_mp_midf)

### ------------------
### Clean All Combined Data

# Further clean the data for modeling

# Filter out plants with no pins (i.e. no seeds found post-dispersal)
kern_w_disp_hw_sv_mp_midf <- kern_w_disp_hw_sv_mp_midf %>% filter(!is.na(Pin))

# Filter out plants with fewer than 3 pins (found seeds) total
Good_IDs <- filter(dplyr::summarize(group_by(kern_w_disp_hw_sv_mp_midf, `Full ID`), npins = n()), npins > 2)$`Full ID`
kern_w_disp_hw_sv_mp_midf <- kern_w_disp_hw_sv_mp_midf %>% filter(`Full ID` %in% Good_IDs)

# Filter out plants with missing maternal length data
kern_w_disp_hw_sv_mp_midf <- filter(kern_w_disp_hw_sv_mp_midf, tot_mat_length_mm > 0.00)

# Change data types 
kern_w_disp_hw_sv_mp_midf$n_inflor <- as.integer(kern_w_disp_hw_sv_mp_midf$n_inflor)
kern_w_disp_hw_sv_mp_midf$SppID <- as.factor(kern_w_disp_hw_sv_mp_midf$SppID)

# clean the column names
colnames(kern_w_disp_hw_sv_mp_midf) <-c("Full_ID", "SppID", "Hub", "Transect", 
                                        "Mom", "Pin", "Distance", "Angle_Comp", 
                                        "Angle_Deg_fromN", "Angle_Rad_fromN", 
                                        "Angle_Deg_fromE", "Angle_Rad_fromE", 
                                        "Notes", "Dist_x", "Dist_y", "mean_Dist_x", 
                                        "mean_Dist_y", "ptcentdist", "Full_ID_dupe", 
                                        "Time", "Day", "mean_speed_1m", "mean_speed_10cm", 
                                        "mean_prop_speed","max_surveg_length_mm", 
                                        "surveg_mass_g", "mean_surveg_length_mm", 
                                        "veg_density_kg_msqr", "mat_mass_mg", 
                                        "main_stem_length_mm", "tot_mat_length_mm", 
                                        "mean_inflor_height_mm", "var_inflor_height_mm", 
                                        "non_na_count", "n_inflor")

kern_w_disp_hw_sv_mp_midf$ptcentdist %>% range

 
# Rescale variables by changing units

# Make sure all distance values have 2 sig figs
kern_w_disp_hw_sv_mp_midf$Distance <- round(kern_w_disp_hw_sv_mp_midf$Distance, digits=2)

# Make sure all dispersion values have 2 sig figs
kern_w_disp_hw_sv_mp_midf$ptcentdist <- round(kern_w_disp_hw_sv_mp_midf$ptcentdist, digits=2)

# Mean proportional wind speed (this is a ratio - scaled to a percentage and decrease sig figs)
kern_w_disp_hw_sv_mp_midf$mean_prop_speed_pcnt <- round(kern_w_disp_hw_sv_mp_midf$mean_prop_speed*100, digits=2)

# Mean surrounding vegetation length (mm -> cm)
kern_w_disp_hw_sv_mp_midf$mean_surveg_length_cm <- round(kern_w_disp_hw_sv_mp_midf$mean_surveg_length_mm/10, digits = 2)

# Density of vegetation surrounding the maternal plant (kg/m^3 -> g/m^3)
kern_w_disp_hw_sv_mp_midf$veg_density_g_msqr <- round(kern_w_disp_hw_sv_mp_midf$veg_density_kg_msqr*1000, digits = 3)

# Main stem length (mm -> cm)
kern_w_disp_hw_sv_mp_midf$main_stem_length_cm <- round(kern_w_disp_hw_sv_mp_midf$main_stem_length_mm/10, digits = 2)

# Total maternal length (mm -> cm)
kern_w_disp_hw_sv_mp_midf$tot_mat_length_cm <- round(kern_w_disp_hw_sv_mp_midf$tot_mat_length_mm/10, digits = 2)

# Mean inflorescence height (mm -> cm)
kern_w_disp_hw_sv_mp_midf$mean_inflor_height_cm <- round(kern_w_disp_hw_sv_mp_midf$mean_inflor_height_mm/10, digits = 2)

# Variance inflorescence height (mm -> cm)
kern_w_disp_hw_sv_mp_midf$var_inflor_height_cm <- round(kern_w_disp_hw_sv_mp_midf$var_inflor_height_mm /100, digits = 2)

# Relative height of inflorescences to surveg
kern_w_disp_hw_sv_mp_midf$rel_inflor_height_svheight <- round(kern_w_disp_hw_sv_mp_midf$mean_inflor_height_cm/ kern_w_disp_hw_sv_mp_midf$mean_surveg_length_cm *100, digits = 2)

# Select columns that are necessary moving forward
kern_w_disp_hw_sv_mp_midf <- dplyr::select(kern_w_disp_hw_sv_mp_midf, "Full_ID", "SppID", "Hub", "Transect", "Mom", "Pin", "Distance",
                                           "ptcentdist", "max_surveg_length_mm", "surveg_mass_g", "mat_mass_mg", "non_na_count", 
                                           "n_inflor", "mean_speed_1m", "mean_speed_10cm", "mean_prop_speed_pcnt", "mean_surveg_length_cm", "veg_density_g_msqr", 
                                           "main_stem_length_cm", "tot_mat_length_cm", "mean_inflor_height_cm","var_inflor_height_cm", "rel_inflor_height_svheight")

### ------------------
### Write Cleaned Data to File

write_csv(kern_w_disp_hw_sv_mp_midf, "./Data/Field Study Cleaned Data - Combined Kernel, Maternal, Sur Veg, and Wind Speed Data.csv")
