#### Wind Tunnel Data Analysis 3: Plot Data and Model Results

# Courtney Van Den Elzen

# Original: 2018 

# Latest Update: January 2023 (syntax update for Jan 2023 version compatibility)

# Description: This script builds scatterplots of each of the predictors against 
#              the distance in the wind tunnel

# WARNING: You must run Wind Tunnel Data Analysis 1 & 2 first to have access to the data

# The seed traits used are:

# 1. Seed mass in mg ("Seed.weight")
# 2. Achene curvature (i.e. angle) measured as max deviance from a straight line ("Achene.twopt.angle.dev.max")
# 3. Ratio of the area of the tip of pappus to the max cross-sectional area of the achene ("pap.tip.cyp.max.area.ratio")
# 4. Ratio of the pappus length to the achene length ("Pap.len.ach.len.ratio")
# 5. Eccentricity (ellipsoidalness) of the pappus tip ("pap.tip.ecc")
# 6. Surface area of the achene ("coneSA")


### ------------------
### Load Packages -------

library(ggplot2)
library(fields)
library(tidyverse)
library(effects)
library(viridis)
library(patchwork)

# Data from Part 2
dist_traits_cal
dist_traits_fre
dist_traits_gla

### ------------------
### Define Custom Functions -------

# This function takes an effect() data frame and calculates the average marginal 
# effect size as a percentage. This is the average because the models are gamma 
# glmms with log links, so the effect size varies across the range of x

PcntEffectVec <-function(effect_df, change_val) {
  unit_eff <- c()
  
  #' Effect range 
  effect_df_rows <- dim(effect_df)[1]
  effect_start <- effect_df$fit[1]
  
  #
  for (i in 1:(dim(effect_df)[1])){
    if (i+change_val <= dim(effect_df)[1]){
      unit_eff[i] <- effect_df[i+change_val, "fit"] - effect_df[i, "fit"]
    }
  }
  pcnt_eff <- unit_eff/effect_start *100
  return(pcnt_eff)
}


UnitEffectVec <-function(effect_df, change_val) {
  unit_eff <- c()
  
  #' 
  for (i in 1:(dim(effect_df)[1])){
    if (i+change_val <= dim(effect_df)[1]){
      unit_eff[i] <- effect_df[i+change_val, "fit"] - effect_df[i, "fit"]
    }
  }
  return(unit_eff)
}


### ------------------
### Data Cleaning and Summaries -------

# Fix the dataframes for plotting and summarizing

# L. californica: select appropriate columns 
dist_traits_cal2 <- dplyr::select(dist_traits_cal, Species.x, Tube.ID, 
                                  Total.Hypoteneuse.Distance, Seed.weight,
                                  Achene.twopt.angle.dev.max, pap.tip.cyp.max.area.ratio, 
                                  Pap.len.ach.len.ratio, pap.tip.ecc, coneSA, 
                                  Helper_Initials)

# Change species coding for plotting
dist_traits_cal2$Species.x <- dplyr::recode(dist_traits_cal2$Species.x, 
                                          `L. californica` = "cal")

# L. fremontii: select appropriate columns 
dist_traits_fre2 <- dplyr::select(dist_traits_fre, Species.x, Tube.ID, Total.Hypoteneuse.Distance,
                           Seed.weight, Achene.twopt.angle.dev.max, pap.tip.cyp.max.area.ratio, 
                           Pap.len.ach.len.ratio, pap.tip.ecc, coneSA, Helper_Initials)

# Change species coding for plotting
dist_traits_fre2$Species.x <- dplyr::recode(dist_traits_fre2$Species.x, `L. fremontii` = "fre")


# L. glaberrima: select appropriate columns 
dist_traits_gla2 <- dplyr::select(dist_traits_gla, Species.x, Tube.ID, Total.Hypoteneuse.Distance,
                           Seed.weight, Achene.twopt.angle.dev.max, pap.tip.cyp.max.area.ratio, 
                           Pap.len.ach.len.ratio, pap.tip.ecc, coneSA, Helper_Initials)

# Change species coding for plotting
dist_traits_gla2$Species.x <- dplyr::recode(dist_traits_gla2$Species.x, `L. glaberrima` = "gla")


# Make a composite dataframe with all species 
dist_traits_cfg <- full_join(dist_traits_cal2, dist_traits_fre2) %>% full_join(dist_traits_gla2) 


### ------------------
### Create Plots  -------

## ------------------
## Seed mass (sm) Plot  -------

# Calculate the marginal effect (slope) and the standard error in the slope for each species
sm_effect_cal <- as.data.frame(effect("Seed.weight", 
                                      mod = calmod_inv, 
                                      xlevels=list(Seed.weight=seq(0.025,0.2,0.025))))

sm_effect_fre <- as.data.frame(effect("Seed.weight", 
                                      mod = fremod_inv, 
                                      xlevels=list(Seed.weight=seq(0.02,0.14,0.02))))

sm_effect_gla <- as.data.frame(effect("Seed.weight", 
                                      mod = glamod_inv, 
                                      xlevels=list(Seed.weight=seq(0.05,0.37,0.02))))

# Calculate the percentage drop in flight distance per 0.1 mg seed mass difference 
# (marginal effect size)
pcnt_eff_sm_cal <- mean(PcntEffectVec(sm_effect_cal, 5))
pcnt_eff_sm_cal

pcnt_eff_sm_fre <- mean(PcntEffectVec(sm_effect_fre, 5))
pcnt_eff_sm_fre

pcnt_eff_sm_gla <- mean(PcntEffectVec(sm_effect_gla, 5))
pcnt_eff_sm_gla

# to get to these "effects" must use mu[distance] = exp(RHS of model)

# Add the species ID as a column
sm_effect_cal$Species.x <- "cal"
sm_effect_fre$Species.x <- "fre"
sm_effect_gla$Species.x <- "gla"

# Create the composite data frame for plotting 
sm_effect <- full_join(sm_effect_cal, sm_effect_fre) %>% full_join(sm_effect_gla) 

# Create the plot
sm_plot <- ggplot(data = sm_effect, aes(Seed.weight, fit, fill = Species.x)) +
  geom_ribbon(aes(ymin=lower, ymax=upper)) + #, alpha= 0.8
  geom_line(colour = "black") +
  geom_point(data = dist_traits_cfg, aes(Seed.weight, 
                                         Total.Hypoteneuse.Distance, 
                                         color = Species.x, 
                                         fill = Species.x)) + #, alpha = 0.4
  scale_colour_manual(values = c("#d95f02", "#21908C", "#440154")) +
  scale_fill_manual(values = c("#d95f02", "#21908C", "#440154")) +
  scale_x_continuous(limits=c(0, 0.4), breaks = seq(0,0.4,0.05)) +
  scale_y_continuous(limits=c(120, 400), breaks = seq(0,400,50)) + 
  labs(
    x = "Seed mass (mg)",
    y = "Distance (cm)"
  ) + 
  guides(fill = "none", colour = "none") +
  theme_classic()

sm_plot

## ------------------
## Achene Curvature  (ca) Plot  -------

# (deviation from straight = 180 degrees)

# Calculate the marginal effect (slope) and the standard error in the slope for each species
ca_effect_cal <- as.data.frame(effect("Achene.twopt.angle.dev.max", 
                                      mod = calmod_inv, 
                                      xlevels=list(Achene.twopt.angle.dev.max=seq(0,23,0.5))))
ca_effect_fre <- as.data.frame(effect("Achene.twopt.angle.dev.max", 
                                      mod = fremod_inv, 
                                      xlevels=list(Achene.twopt.angle.dev.max=seq(0,30.5,0.5))))
ca_effect_gla <- as.data.frame(effect("Achene.twopt.angle.dev.max", 
                                      mod = glamod_inv, 
                                      xlevels=list(Achene.twopt.angle.dev.max=seq(0,17,0.5))))

#' Add the species ID as a column
ca_effect_cal$Species.x <- "cal"
ca_effect_fre$Species.x <- "fre"
ca_effect_gla$Species.x <- "gla"

#' Calculate the percentage drop in flight distance per 1 unit difference (marginal effect size)
pcnt_eff_ca_cal <- mean(PcntEffectVec(ca_effect_cal, 2))
pcnt_eff_ca_cal

pcnt_eff_ca_fre <- mean(PcntEffectVec(ca_effect_fre, 2))
pcnt_eff_ca_fre

pcnt_eff_ca_gla <- mean(PcntEffectVec(ca_effect_gla, 2))
pcnt_eff_ca_gla

#' Create the composite data frame for plotting
ca_effect <- full_join(ca_effect_cal, ca_effect_fre) %>% full_join(ca_effect_gla)

#' Create the plot
ca_plot <- ggplot(data = ca_effect, aes(Achene.twopt.angle.dev.max, fit, fill = Species.x)) +
  geom_ribbon(aes(ymin=lower, ymax=upper)) + #, alpha = 0.8
  geom_line(colour = "black")+
  geom_point(data = dist_traits_cfg, aes(Achene.twopt.angle.dev.max, Total.Hypoteneuse.Distance, color = Species.x, fill = Species.x))+  #, alpha = 0.4
  scale_colour_manual(values = c("#d95f02", "#21908C", "#440154")) +
  scale_fill_manual(values = c("#d95f02", "#21908C", "#440154")) +
  scale_x_continuous(limits=c(0, 31), breaks = seq(0,31,5)) +
  scale_y_continuous(limits=c(120, 400), breaks = seq(0,400,50)) + 
  labs(
    x = "Angle of the cypsela (°)",
    y = ""
  ) + 
  theme(axis.title.y = element_blank()) +
  theme_classic() +
  guides(fill = "none", colour = "none") 

ca_plot


## ------------------
## Pappus Area to Cypsela Area Ratio (ar) Plot  -------

# Calculate the marginal effect (slope) and the standard error in the slope for each species
ar_effect_cal <- as.data.frame(effect("pap.tip.cyp.max.area.ratio", 
                                      mod = calmod_inv, 
                                      xlevels=list(pap.tip.cyp.max.area.ratio=seq(0.5,30,0.5))))
ar_effect_fre <- as.data.frame(effect("pap.tip.cyp.max.area.ratio", 
                                      mod = fremod_inv, 
                                      xlevels=list(pap.tip.cyp.max.area.ratio=seq(1,19,0.5))))
ar_effect_gla <- as.data.frame(effect("pap.tip.cyp.max.area.ratio", 
                                      mod = glamod_inv, 
                                      xlevels=list(pap.tip.cyp.max.area.ratio=seq(0.5,5,0.05))))

# Add the species ID as a column
ar_effect_cal$Species.x <- "cal"
ar_effect_fre$Species.x <- "fre"
ar_effect_gla$Species.x <- "gla"

# Create the composite data frame for plotting
ar_effect <- full_join(ar_effect_cal, ar_effect_fre) %>% full_join(ar_effect_gla)

# Calculate the percentage drop in flight distance per 1 unit difference (marginal effect size)
pcnt_eff_ar_cal <- mean(PcntEffectVec(ar_effect_cal, 2))
pcnt_eff_ar_cal

pcnt_eff_ar_fre <- mean(PcntEffectVec(ar_effect_fre, 2))
pcnt_eff_ar_fre

pcnt_eff_ar_gla <- mean(PcntEffectVec(ar_effect_gla, 20))
pcnt_eff_ar_gla

# Create the plot
ar_plot <- ggplot(data = ar_effect, aes(pap.tip.cyp.max.area.ratio, fit, fill = Species.x)) +
  geom_ribbon(aes(ymin=lower, ymax=upper)) + #, alpha = 0.8
  geom_line(colour="black")+
  geom_point(data = dist_traits_cfg, aes(pap.tip.cyp.max.area.ratio, Total.Hypoteneuse.Distance, color = Species.x, fill = Species.x))+ #, alpha = 0.4
  scale_colour_manual(values = c("#d95f02", "#21908C", "#440154"), labels = c("L.californica","L. fremontii","L. glaberrima")) +
  scale_fill_manual(values = c("#d95f02", "#21908C", "#440154"), labels = c("L.californica","L. fremontii","L. glaberrima")) +
  scale_x_continuous(limits=c(0, 30), breaks = seq(0,31,5)) +
  scale_y_continuous(limits=c(120, 400), breaks = seq(0,400,50)) + 
  labs(
    x = "Pappus tip area (mm^2)/cypsela cs area (mm^2)",
    y = ""
  ) + 
  theme_classic() + 
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(face = "italic"))

ar_plot


## ------------------
## Pappus Length to Cypsela Length Ratio (lr) Plot  -------

# Calculate the marginal effect (slope) and the standard error in the slope for each species
lr_effect_cal <- as.data.frame(effect("Pap.len.ach.len.ratio", 
                                      mod = calmod_inv, 
                                      xlevels=list(Pap.len.ach.len.ratio=seq(0,1.97,0.02))))
lr_effect_fre <- as.data.frame(effect("Pap.len.ach.len.ratio", 
                                      mod = fremod_inv, 
                                      xlevels=list(Pap.len.ach.len.ratio=seq(1.06,2.10,0.02))))
lr_effect_gla <- as.data.frame(effect("Pap.len.ach.len.ratio", 
                                      mod = glamod_inv, 
                                      xlevels=list(Pap.len.ach.len.ratio=seq(0.27,0.67,0.02))))

#' Calculate the percentage drop in flight distance per 1 unit difference (marginal effect size)
pcnt_eff_lr_cal <- mean(PcntEffectVec(lr_effect_cal, 6))
pcnt_eff_lr_cal

pcnt_eff_lr_fre <- mean(PcntEffectVec(lr_effect_fre, 6))
pcnt_eff_lr_fre

pcnt_eff_lr_gla <- mean(PcntEffectVec(lr_effect_gla, 6))
pcnt_eff_lr_gla

#' Add the species ID as a column
lr_effect_cal$Species.x <- "cal"
lr_effect_fre$Species.x <- "fre"
lr_effect_gla$Species.x <- "gla"

#' Create the composite data frame for plotting
lr_effect <- full_join(lr_effect_cal, lr_effect_fre) %>% full_join(lr_effect_gla)

#' Create the plot
lr_plot <- ggplot(data = lr_effect, aes(Pap.len.ach.len.ratio, fit, fill = Species.x)) +
  geom_ribbon(aes(ymin=lower, ymax=upper)) + #, alpha = 0.8
  geom_line(colour = "black")+
  geom_point(data = dist_traits_cfg, aes(Pap.len.ach.len.ratio, Total.Hypoteneuse.Distance, color = Species.x, fill = Species.x))+ #, alpha = 0.4
  scale_colour_manual(values = c("#d95f02", "#21908C", "#440154")) +
  scale_fill_manual(values = c("#d95f02", "#21908C", "#440154")) +
  scale_y_continuous(limits=c(120, 400), breaks = seq(0,400,50)) + 
  scale_x_continuous(limits=c(0, 2.25), breaks = seq(0,2.25,0.5)) +
  labs(
    x = "Pappus length (mm)/cypsela length (mm)",
    y = "Distance (cm)"
  ) + 
  guides(fill = "none", colour = "none")+
  theme_classic()

lr_plot


## ------------------
## Pappus Tip Eccentricity (pe) Plot  -------

# Calculate the marginal effect (slope) and the standard error in the slope for each species
pe_effect_cal <- as.data.frame(effect("pap.tip.ecc", 
                                      mod = calmod_inv, 
                                      xlevels=list(pap.tip.ecc=seq(0,0.99,0.01))))
pe_effect_fre <- as.data.frame(effect("pap.tip.ecc", 
                                      mod = fremod_inv, 
                                      xlevels=list(pap.tip.ecc=seq(0,0.92,0.01))))
pe_effect_gla <- as.data.frame(effect("pap.tip.ecc", 
                                      mod = glamod_inv, 
                                      xlevels=list(pap.tip.ecc=seq(0,0.88,0.01))))

# Add the species ID as a column
pe_effect_cal$Species.x <- "cal"
pe_effect_fre$Species.x <- "fre"
pe_effect_gla$Species.x <- "gla"

# Create the composite data frame for plotting
pe_effect <- full_join(pe_effect_cal, pe_effect_fre) %>% full_join(pe_effect_gla)

#' Create the plot
pe_plot <- ggplot(data = pe_effect, aes(pap.tip.ecc, fit, fill = Species.x)) +
  geom_ribbon(aes(ymin=lower, ymax=upper)) + #, alpha = 0.8
  geom_line(colour = "black")+
  geom_point(data = dist_traits_cfg, aes(pap.tip.ecc, Total.Hypoteneuse.Distance, color = Species.x, fill = Species.x))+ #, alpha = 0.4
  scale_colour_manual(values = c("#d95f02", "#21908C", "#440154")) +
  scale_fill_manual(values = c("#d95f02", "#21908C", "#440154")) +
  scale_y_continuous(limits=c(120, 400), breaks = seq(0,400,50)) + 
  labs(
    x = "Pappus tip eccentricity",
    y = ""
  ) + 
  theme_classic() +
  guides(fill = "none", colour = "none")

pe_plot


## ------------------
## Surface Area of the Cypsela (sa) Plot  -------

# Calculate the marginal effect (slope) and the standard error in the slope for each species
sa_effect_cal <- as.data.frame(effect("coneSA", 
                                      mod = calmod_inv, 
                                      xlevels=list(coneSA=seq(0.4,1.16,0.01))))
sa_effect_fre <- as.data.frame(effect("coneSA", 
                                      mod = fremod_inv, 
                                      xlevels=list(coneSA=seq(0.3,0.75,0.01))))
sa_effect_gla <- as.data.frame(effect("coneSA", 
                                      mod = glamod_inv, 
                                      xlevels=list(coneSA=seq(0.6,1.82,0.01))))

# Add the species ID as a column
sa_effect_cal$Species.x <- "cal"
sa_effect_fre$Species.x <- "fre"
sa_effect_gla$Species.x <- "gla"

# Create the composite data frame for plotting
sa_effect <- full_join(sa_effect_cal, sa_effect_fre) %>% full_join(sa_effect_gla)

# Create the plot
sa_plot <- ggplot(data = sa_effect, aes(coneSA, fit, fill = Species.x)) +
  geom_ribbon(aes(ymin=lower, ymax=upper)) + #, alpha = 0.8
  geom_line(colour = "black")+
  geom_point(data = dist_traits_cfg, aes(coneSA, Total.Hypoteneuse.Distance, color = Species.x, fill = Species.x))+ # alpha = 0.4
  scale_colour_manual(values = c("#d95f02", "#21908C", "#440154")) +
  scale_fill_manual(values = c("#d95f02", "#21908C", "#440154")) +
  scale_y_continuous(limits=c(120, 400), breaks = seq(0,400,50)) + 
  scale_x_continuous(limits=c(0.25,2.0), breaks = seq(0.25,2.0,0.25)) +
  labs(
    x = "Cypsela surface area (mm^2)",
    y = ""
  ) + 
  theme_classic() +
  guides(fill = "none", colour = "none")

sa_plot


## ------------------
## Combined Panel Plot  -------

# Use the patchwork package to make this plot

layout <- c(
  area(1, 2, 1, 5),
  area(2, 1, 3, 2),
  area(2, 3, 3, 4),
  area(2, 5, 3, 6),
  area(4, 1, 5, 2),
  area(4, 3, 5, 4),
  area(4, 5, 5, 6))

composite_plot <- 
  guide_area() + 
  lr_plot + ar_plot + pe_plot + sm_plot + ca_plot + sa_plot + 
  plot_layout(design = layout, guides = "collect")

composite_plot

