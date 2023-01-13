#### 2019 Field Dispersal Data Analysis 7: Inter-seed Spread Scatterplots 

# Courtney Van Den Elzen

# Initiated: February 2020
# Latest Update: January 2023 (syntax update for Jan 2023 version compatibility)

# Description: Building the inter-seed spread scatterplots (Figure 4, S4C,D)

# WARNING: You must first run Field Study Data Analysis 5 - Modeling.R so the 
#          plots have access to the models


### ------------------
### Load Packages -------
library(ggplot2)
library(fields)
library(tidyverse)
library(effects)
library(viridis)
library(patchwork)

### ------------------
### Clean Data ------

# change data types for plotting function compatibility
kern_hw_sv_mp_midf <- kern_hw_sv_mp_midf %>% mutate(
  SppID = factor(SppID), 
  Full_ID = factor(Full_ID), 
  Hub = factor(Hub),
  Transect = factor(Transect),
  Mom = factor(Mom),
  Pin = factor(Pin)
)

kern_hw_sv_mp_midf$SppID <- as.factor(kern_hw_sv_mp_midf$SppID)
kern_hw_sv_mp_midf_frebase <- within(kern_hw_sv_mp_midf, SppID <- relevel(SppID, ref = "fre"))

# Split by species
modsel_cal <- filter(kern_hw_sv_mp_midf, SppID == "cal")
modsel_fre <- filter(kern_hw_sv_mp_midf, SppID == "fre")
modsel_gla <- filter(kern_hw_sv_mp_midf, SppID == "gla")

# Fix the dataframes for plotting and summarizing

modsel_cal2 <- dplyr::select(modsel_cal, Full_ID, SppID, Distance, 
                      ptcentdist, n_inflor, mean_prop_speed_pcnt, 
                      mean_surveg_length_cm, veg_density_g_msqr, 
                      main_stem_length_cm, tot_mat_length_cm, 
                      mean_inflor_height_cm, var_inflor_height_cm)

modsel_fre2 <- dplyr::select(modsel_fre, Full_ID, SppID, Distance, 
                      ptcentdist, n_inflor, mean_prop_speed_pcnt, 
                      mean_surveg_length_cm, veg_density_g_msqr, 
                      main_stem_length_cm, tot_mat_length_cm, 
                      mean_inflor_height_cm, var_inflor_height_cm)

modsel_gla2 <- dplyr::select(modsel_gla, Full_ID, SppID, Distance, 
                      ptcentdist, n_inflor, mean_prop_speed_pcnt, 
                      mean_surveg_length_cm, veg_density_g_msqr, 
                      main_stem_length_cm, tot_mat_length_cm, 
                      mean_inflor_height_cm, var_inflor_height_cm)

# Make some composite dataframes

modsel_c_mean <- group_by(modsel_cal2, Full_ID, SppID, n_inflor, mean_prop_speed_pcnt, 
                          mean_surveg_length_cm, veg_density_g_msqr, main_stem_length_cm, 
                          tot_mat_length_cm, mean_inflor_height_cm, var_inflor_height_cm) %>%
  dplyr::summarize(mean_dist = mean(Distance), 
                   mean_disp = mean(ptcentdist), 
                   se_dist = se(Distance), 
                   se_disp = se(ptcentdist),
                   min_dist = min(Distance),
                   max_dist = max(Distance),
                   min_disp = min(ptcentdist),
                   max_disp = max(ptcentdist))

modsel_f_mean <- group_by(modsel_fre2, Full_ID, SppID, n_inflor, mean_prop_speed_pcnt, 
                          mean_surveg_length_cm, veg_density_g_msqr, main_stem_length_cm, 
                          tot_mat_length_cm, mean_inflor_height_cm, var_inflor_height_cm) %>%
  dplyr::summarize(mean_dist = mean(Distance), 
                   mean_disp = mean(ptcentdist), 
                   se_dist = se(Distance), 
                   se_disp = se(ptcentdist),
                   min_dist = min(Distance),
                   max_dist = max(Distance),
                   min_disp = min(ptcentdist),
                   max_disp = max(ptcentdist))

modsel_g_mean <- group_by(modsel_gla2, Full_ID, SppID, n_inflor, mean_prop_speed_pcnt, 
                          mean_surveg_length_cm, veg_density_g_msqr, main_stem_length_cm, 
                          tot_mat_length_cm, mean_inflor_height_cm, var_inflor_height_cm) %>%
  dplyr::summarize(mean_dist = mean(Distance), 
                   mean_disp = mean(ptcentdist), 
                   se_dist = se(Distance), 
                   se_disp = se(ptcentdist),
                   min_dist = min(Distance),
                   max_dist = max(Distance),
                   min_disp = min(ptcentdist),
                   max_disp = max(ptcentdist))

modsel_cg_mean <- full_join(modsel_cal2, modsel_gla2) %>% 
  group_by(Full_ID, SppID, n_inflor, mean_prop_speed_pcnt, 
           mean_surveg_length_cm, veg_density_g_msqr, 
           main_stem_length_cm, tot_mat_length_cm, 
           mean_inflor_height_cm, var_inflor_height_cm) %>%
  dplyr::summarize(mean_dist = mean(Distance),
                   mean_disp = mean(ptcentdist), 
                   se_dist = se(Distance), 
                   se_disp = se(ptcentdist),
                   min_dist = min(Distance),
                   max_dist = max(Distance),
                   min_disp = min(ptcentdist),
                   max_disp = max(ptcentdist))

modsel_cfg_mean <- full_join(modsel_cal2, modsel_fre2) %>% 
  full_join(modsel_gla2) %>%
  group_by(Full_ID, SppID, n_inflor, mean_prop_speed_pcnt, 
           mean_surveg_length_cm, veg_density_g_msqr, 
           main_stem_length_cm, tot_mat_length_cm, 
           mean_inflor_height_cm, var_inflor_height_cm) %>%
  dplyr::summarize(mean_dist = mean(Distance),
                   mean_disp = mean(ptcentdist), 
                   se_dist = se(Distance), 
                   se_disp = se(ptcentdist),
                   min_dist = min(Distance),
                   max_dist = max(Distance),
                   min_disp = min(ptcentdist),
                   max_disp = max(ptcentdist))


## ------
## Mean Inflorescence Height Plot -----

# Calculate the marginal effect (slope) and the standard error in the slope for each species
ih_effect_cal_disp <- as.data.frame(effect("mean_inflor_height_cm", mod = cal_id_disp, xlevels=list(mean_inflor_height_cm=seq(7,23,1))))
ih_effect_fre_disp <- as.data.frame(effect("mean_inflor_height_cm", mod = fre_id_disp, xlevels=list(mean_inflor_height_cm=seq(12,21,1))))
ih_effect_gla_disp <- as.data.frame(effect("mean_inflor_height_cm", mod = gla_id_disp, xlevels=list(mean_inflor_height_cm=seq(17,35,1))))

# Add the species ID as a column
ih_effect_cal_disp$SppID <- "cal"
ih_effect_fre_disp$SppID <- "fre"
ih_effect_gla_disp$SppID <- "gla"

# Create the composite data frame for plotting
ih_effect_disp <- full_join(ih_effect_cal_disp, ih_effect_fre_disp) %>% full_join(ih_effect_gla_disp)

# Create the plot
ih_plot_disp <- ggplot(data = ih_effect_disp, aes(mean_inflor_height_cm, fit, fill = SppID)) +
  geom_ribbon(aes(ymin=lower, ymax=upper)) + #, alpha = 0.8
  geom_line(colour="black")+
  geom_point(data = modsel_cfg_mean, aes(mean_inflor_height_cm, mean_disp, color = SppID, fill = SppID), shape = 20, size = 3)+  #, alpha = 0.4
  geom_point(data = modsel_cal2, aes(mean_inflor_height_cm, ptcentdist), color = "#d95f02", shape=20, size = 0.3) + #, alpha = 0.4
  geom_point(data = modsel_fre2, aes(mean_inflor_height_cm, ptcentdist), color = "#21908C", shape=20, size = 0.3) + #, alpha = 0.4
  geom_point(data = modsel_gla2, aes(mean_inflor_height_cm, ptcentdist), color = "#440154", shape=20, size = 0.3) + #, alpha = 0.4
  scale_colour_manual(values = c("#d95f02", "#21908C", "#440154"), labels = c("L.californica","L. fremontii","L. glaberrima")) +
  scale_fill_manual(values = c("#d95f02", "#21908C", "#440154"), labels = c("L.californica","L. fremontii","L. glaberrima")) +
  labs(
    x = "Mean infloresence height (cm)",
    y = ""
  ) + 
  scale_y_continuous(limits=c(-1, 25), breaks = seq(0,25,5)) +
  theme_classic() + 
  theme(axis.title.y = element_blank(),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(face = "italic"))

ih_plot_disp

## ------
## Number of Inflorescences Plot -----
 
# Calculate the marginal effect (slope) and the standard error in the slope for each species
ni_effect_cal_disp <- as.data.frame(effect("n_inflor", mod = cal_id_disp, xlevels=list(n_inflor=seq(1,3,1))))
ni_effect_gla_disp <- as.data.frame(effect("n_inflor", mod = gla_id_disp, xlevels=list(n_inflor=seq(1,7,1))))

# Add the species ID as a column
ni_effect_cal_disp$SppID <- "cal"
ni_effect_gla_disp$SppID <- "gla"

# Create the composite data frame for plotting
ni_effect_disp <- full_join(ni_effect_cal_disp, ni_effect_gla_disp)

# Create the plot
ni_plot_disp <- ggplot(data = ni_effect_disp, aes(n_inflor, fit, fill = SppID)) +
  geom_ribbon(aes(ymin=lower, ymax=upper)) + #, alpha = 0.8
  geom_line(colour = "black")+
  geom_point(data = modsel_cg_mean, aes(n_inflor, mean_disp, color = SppID, fill = SppID), shape = 20, size = 3)+ #, alpha = 0.4
  geom_point(data = modsel_cal2, aes(n_inflor, ptcentdist), color = "#d95f02", shape=20, size = 0.3) + #, alpha = 0.4
  geom_point(data = modsel_gla2, aes(n_inflor, ptcentdist), color = "#440154", shape=20, size = 0.3) + #, alpha = 0.4
  scale_colour_manual(values = c("#d95f02", "#440154")) +
  scale_fill_manual(values = c("#d95f02", "#440154")) + 
  labs(
    x = "Number of infloresences",
    y = "Dispersion (cm)"
  ) + 
  scale_y_continuous(limits=c(-1, 25), breaks = seq(0,25,5)) +
  guides(fill = "none", colour = "none") +
  theme_classic()

ni_plot_disp

## ------
## Variance in Inflorescence Height Plot -----

#' Calculate the marginal effect (slope) and the standard error in the slope for each species
vih_effect_gla_disp <- as.data.frame(effect("var_inflor_height_cm", mod = gla_id_disp, xlevels=list(var_inflor_height_cm=seq(0,20,2))))
vih_effect_gla_disp

#' Add the species ID as a column
vih_effect_gla_disp$SppID <- "gla"

#' Create the composite data frame for plotting
vih_effect_disp <- vih_effect_gla_disp

#' Create the plot
vih_plot_disp <- ggplot(data = vih_effect_disp, aes(var_inflor_height_cm, fit, fill = SppID)) +
  geom_ribbon(aes(ymin=lower, ymax=upper)) + #, alpha = 0.8
  geom_line(colour = "black")+
  geom_point(data = modsel_g_mean, aes(var_inflor_height_cm, mean_disp, color = SppID, fill = SppID), shape = 20, size = 3)+ #, alpha = 0.4
  geom_point(data = modsel_gla2, aes(var_inflor_height_cm, ptcentdist), color = "#440154", shape=20, size = 0.3) + #, alpha = 0.4
  scale_colour_manual(values = c("#440154")) +
  scale_fill_manual(values = c("#440154")) + 
  labs(
    x = "Variance in inflorescence height (cm^2)",
    y = ""
  ) + 
  theme_classic() +
  theme(axis.title.y = element_blank()) +
  scale_y_continuous(limits=c(-1, 25), breaks = seq(0,25,5)) +
  guides(fill = "none", colour = "none")

vih_plot_disp

## ------
## Wind Exposure Covariate Plot -----

# Calculate the marginal effect (slope) and the standard error in the slope for each species
ws_effect_cal_disp <- as.data.frame(effect("mean_prop_speed_pcnt", mod = cal_id_disp, xlevels=list(mean_prop_speed_pcnt=seq(7,46,2))))
ws_effect_fre_disp <- as.data.frame(effect("mean_prop_speed_pcnt", mod = fre_id_disp, xlevels=list(mean_prop_speed_pcnt=seq(9,52,2))))
ws_effect_gla_disp <- as.data.frame(effect("mean_prop_speed_pcnt", mod = gla_id_disp, xlevels=list(mean_prop_speed_pcnt=seq(0,64,2))))

# Add the species ID as a column
ws_effect_cal_disp$SppID <- "cal"
ws_effect_fre_disp$SppID <- "fre"
ws_effect_gla_disp$SppID <- "gla"

# Create the composite data frame for plotting 
ws_effect_disp <- full_join(ws_effect_cal_disp, ws_effect_fre_disp) %>% full_join(ws_effect_gla_disp)

# Create the plot
ws_plot_disp <- ggplot(data = ws_effect_disp, aes(mean_prop_speed_pcnt, fit, fill = SppID)) +
  geom_ribbon(aes(ymin=lower, ymax=upper)) + #, alpha = 0.8
  geom_line(colour = "black") +
  geom_point(data = modsel_cfg_mean, aes(mean_prop_speed_pcnt, mean_disp, color = SppID, fill = SppID), shape = 20, size = 3) + #, alpha = 0.4
  geom_point(data = modsel_cal2, aes(mean_prop_speed_pcnt, ptcentdist), color = "#d95f02", shape=20, size = 0.3) + #, alpha = 0.4
  geom_point(data = modsel_fre2, aes(mean_prop_speed_pcnt, ptcentdist), color = "#21908C", shape=20, size = 0.3) + #, alpha = 0.4
  geom_point(data = modsel_gla2, aes(mean_prop_speed_pcnt, ptcentdist), color = "#440154", shape=20, size = 0.3) + #, alpha = 0.4
  scale_colour_manual(values = c("#d95f02", "#21908C", "#440154")) +
  scale_fill_manual(values = c("#d95f02", "#21908C", "#440154")) + 
  labs(
    x = "Mean proportional wind speed (%)",
    y = "Inter-seed spread (cm)"
  ) + 
  scale_y_continuous(limits=c(-1, 25), breaks = seq(0,25,5)) +
  guides(fill = "none", colour = "none") +
  theme_classic()

ws_plot_disp

## ------
## Surrounding Vegetation Length Covariate Plot -----

# Calculate the marginal effect (slope) and the standard error in the slope for each species
sv_effect_cal_disp <- as.data.frame(effect("mean_surveg_length_cm", mod = cal_id_disp, xlevels=list(mean_surveg_length_cm=seq(11,30,1))))
sv_effect_fre_disp <- as.data.frame(effect("mean_surveg_length_cm", mod = fre_id_disp, xlevels=list(mean_surveg_length_cm=seq(14,31,1))))
sv_effect_gla_disp <- as.data.frame(effect("mean_surveg_length_cm", mod = gla_id_disp, xlevels=list(mean_surveg_length_cm=seq(16,39,1))))

#' Add the species ID as a column
sv_effect_cal_disp$SppID <- "cal"
sv_effect_fre_disp$SppID <- "fre"
sv_effect_gla_disp$SppID <- "gla"

#' Create the composite data frame for plotting
sv_effect_disp <- full_join(sv_effect_cal_disp, sv_effect_fre_disp) %>% full_join(sv_effect_gla_disp)

#' Create the plot
sv_plot_disp <- ggplot(data = sv_effect_disp, aes(mean_surveg_length_cm, fit, fill = SppID)) +
  geom_ribbon(aes(ymin=lower, ymax=upper)) +  #, alpha = 0.8
  geom_line(colour = "black")+
  geom_point(data = modsel_cfg_mean, aes(mean_surveg_length_cm, mean_disp, color = SppID, fill = SppID), shape = 20, size = 3)+  #, alpha = 0.4
  geom_point(data = modsel_cal2, aes(mean_surveg_length_cm, ptcentdist), color = "#d95f02", shape=20, size = 0.3) + #, alpha = 0.4
  geom_point(data = modsel_fre2, aes(mean_surveg_length_cm, ptcentdist), color = "#21908C", shape=20, size = 0.3) + #, alpha = 0.4
  geom_point(data = modsel_gla2, aes(mean_surveg_length_cm, ptcentdist), color = "#440154", shape=20, size = 0.3) + #, alpha = 0.4
  scale_colour_manual(values = c("#d95f02", "#21908C", "#440154")) +
  scale_fill_manual(values = c("#d95f02", "#21908C", "#440154")) +
  labs(
    x = "Mean surrounding vegetation height (cm)",
    y = ""
  ) + 
  theme(axis.title.y = element_blank()) +
  scale_y_continuous(limits=c(-5, 25), breaks = seq(0,25,5)) +
  theme_classic() +
  guides(fill = "none", colour = "none") 

sv_plot_disp

## ------
## Combine All Plots Into a Panel Plot -----

# Use the patchwork package to combine the plots into 1. 
layout <- c(
  area(1, 2, 1, 5),
  area(2, 1, 3, 2),
  area(2, 3, 3, 4),
  area(2, 5, 3, 6),
  area(4, 2, 5, 3),
  area(4, 4, 5, 5))

composite_plot_disp <- 
  guide_area() + 
  ni_plot_disp + 
  ih_plot_disp + 
  vih_plot_disp + 
  ws_plot_disp + 
  sv_plot_disp + 
  plot_layout(design = layout, guides = "collect")

composite_plot_disp
