#### 2019 Field Dispersal Data Analysis 6: Density Plots by Species

# Courtney Van Den Elzen

# Initiated: February 2020
# Latest Update: January 2023 (syntax update for Jan 2023 version compatibility)

# Description: Building the density plots (Figure 2)


### ------------------
### Load Packages -------

library(ggridges)
library(ggplot2)
library(viridis)
library(tidyverse)

### ------------------
### Read In Data -------

kern_hw_sv_mp_midf <- read_csv("./Data/Field Study Cleaned Data - Combined Kernel, Maternal, Sur Veg, and Wind Speed Data.csv")

### ------------------
### Clean Data 

# change the data types for compatibility with plot functions
kern_hw_sv_mp_midf <- kern_hw_sv_mp_midf %>% mutate(
  SppID = factor(SppID), 
  Full_ID = factor(Full_ID), 
  Hub = factor(Hub),
  Transect = factor(Transect),
  Mom = factor(Mom),
  Pin = factor(Pin)
)

# Break the data down by species
modsel_cal <- filter(kern_hw_sv_mp_midf, SppID == "cal")
modsel_fre <- filter(kern_hw_sv_mp_midf, SppID == "fre")
modsel_gla <- filter(kern_hw_sv_mp_midf, SppID == "gla")

# Relevel the species so they're in the desired order
kern_hw_sv_mp_midf$SppID <- relevel(kern_hw_sv_mp_midf$SppID, "fre")
kern_hw_sv_mp_midf$SppID <- relevel(kern_hw_sv_mp_midf$SppID, "gla")
levels(kern_hw_sv_mp_midf$SppID) <- c("L. glaberrima", "L. fremontii", "L. californica")

# Make species-specific dataframes
kern_hw_sv_mp_midf_cal <- kern_hw_sv_mp_midf %>% dplyr::filter(SppID == "L. californica")
kern_hw_sv_mp_midf_fre <- kern_hw_sv_mp_midf %>% dplyr::filter(SppID == "L. fremontii")
kern_hw_sv_mp_midf_gla <- kern_hw_sv_mp_midf %>% dplyr::filter(SppID == "L. glaberrima")

# Data summaries: distance (mean and se)
meandf_dist <- 
  kern_hw_sv_mp_midf %>% 
  group_by(Full_ID, SppID) %>% 
  dplyr::summarize(mean_dist = mean(Distance), n_dist = n()) %>% 
  group_by(SppID) %>%  
  dplyr::summarize(mean_mean_dist = mean(mean_dist), 
                   # standard error
                   se_mean_dist = sd(mean_dist)/sqrt(length(mean_dist)))

# Data summaries: inter-seed spread (mean and se)
meandf_disp <- 
  kern_hw_sv_mp_midf %>% 
  group_by(Full_ID, SppID) %>% 
  dplyr::summarize(mean_disp = mean(ptcentdist)) %>% 
  group_by(SppID)%>% 
  dplyr::summarize(mean_mean_disp = mean(mean_disp), 
                   se_mean_disp = sd(mean_disp)/sqrt(length(mean_disp)))


### ------------------
### Figure 2: Density Plots -------

# using `ggridges` to create plot of the distance from the base of the maternal plant
dist_panel_density <- ggplot(kern_hw_sv_mp_midf, aes(x = Distance, y = SppID, fill = SppID)) +
  geom_density_ridges(
    jittered_points = FALSE,
    quantile_lines = FALSE,
    scale = 0.95,
    alpha = 0.7) + # adding annotations
  geom_vline(xintercept=c(meandf_dist$mean_mean_dist), 
             colour = viridis(3), 
             linetype = "dashed",
             size = 0.8) + # L. cal mean
  geom_point(aes(x=meandf_dist$mean_mean_dist[3], y=3),shape=21, size = 2, fill = rev(viridis(3))[1]) + # L. cal mean
  geom_point(aes(x=meandf_dist$mean_mean_dist[2], y=2),shape=21, size = 2, fill = rev(viridis(3))[2]) + # L. fre mean
  geom_point(aes(x=meandf_dist$mean_mean_dist[1], y=1.01),shape=21, size = 2, fill = rev(viridis(3))[3]) + # L. gla mean
  scale_x_continuous(expand = c(0.005, 0), breaks = c(0,5,10,15,20,25,30,35)) +
  labs(
    x = "Distance from maternal plant (cm)",
    y = "Density"
  ) + 
  theme_classic() +
  scale_y_discrete(expand = c(0,0)) +
  theme(legend.position = "none", 
        axis.text.y = element_text(face = "italic"),
        text=element_text(size=20)) +
  scale_fill_viridis(discrete=TRUE, direction = 1)


# using `ggridges` to create plot of the inter-seed spread
disp_panel_density <- ggplot(kern_hw_sv_mp_midf, aes(x = ptcentdist, y = SppID, fill = SppID)) +
  geom_density_ridges(
    jittered_points = FALSE,
    quantile_lines = FALSE,
    scale = 0.95,
    alpha = 0.7) + # adding annotation
  geom_vline(xintercept=meandf_disp$mean_mean_disp,
             colour = viridis(3),
             linetype = "dashed",
             size=0.8) + # L. cal mean
  geom_point(aes(x=meandf_disp$mean_mean_disp[3], y=3),shape=21, size = 2, fill = rev(viridis(3))[1]) + # L. cal mean
  geom_point(aes(x=meandf_disp$mean_mean_disp[2], y=2),shape=21, size = 2, fill = rev(viridis(3))[2]) + # L. fre mean
  geom_point(aes(x=meandf_disp$mean_mean_disp[1], y=1.01),shape=21, size = 2, fill = rev(viridis(3))[3]) + # L. gla mean
  scale_y_discrete(expand = c(0,0)) + 
  scale_x_continuous(expand = c(0.005, 0), breaks = c(0,5,10,15,20,25,30)) +
  labs(
    x = "Inter-seed Spread (cm)",
    y = "Density"
  ) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.y = element_text(face = "italic"),
        text=element_text(size=20)) +
  scale_fill_viridis(direction=1, discrete=TRUE)



