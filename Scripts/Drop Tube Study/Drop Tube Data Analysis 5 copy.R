#' # Drop Tube Data Analysis Part 5
#' 
#' ### Create the panel plot of the scatter plots of each predictor vs the response variable
#' 
#' Courtney Van Den Elzen 
#'
#' 2018-2019 Study
#' 

#' #### Libraries
library(lme4)
library(tidyverse)
library(nlme)
library(cubature)
library(glue)
library(car)
library(fields)
library(effects)
library(viridis)
library(patchwork)

#' Here is the data
dt_traits_merge_cal
dt_traits_merge_fre
dt_traits_merge_gla

#' Custom functions
#' 
#' This function takes an effect() data frame and calculates the average marginal effect size as a percentage.
#' This is the average because the models are gamma glmms with log links, so the effect size varies across the range of x
PcntEffectVec <-function(effect_df, change_val) {
  unit_eff <- c()
  
  #' Effect range 
  effect_df_rows <- dim(effect_df)[1]
  effect_start <- effect_df$fit[1]
  
  #' 
  for (i in 1:(dim(effect_df)[1])){
    if (i+change_val <= dim(effect_df)[1]){
      unit_eff[i] <- effect_df[i+change_val, "fit"] - effect_df[i, "fit"]
    }
  }
  pcnt_eff <- unit_eff/effect_start *100
  return(pcnt_eff)
}

#' #### Fix the dataframes for plotting and summarizing
dt_traits_merge_cal2 <- dplyr::select(dt_traits_merge_cal, Velocity_Vals_m_sec, Seed_Mass_mg, 
                                 Seed_Type, Achene.twopt.angle.dev.max, pap.tip.cyp.max.area.ratio,
                                 Pap.len.ach.len.ratio, pap.tip.ecc, coneSA, Plant_ID, Species)

dt_traits_merge_cal2$Species <- dplyr::recode(dt_traits_merge_cal2$Species, `L. californica` = "cal")

dt_traits_merge_fre2 <- dplyr::select(dt_traits_merge_fre, Velocity_Vals_m_sec, Seed_Mass_mg, 
                               Seed_Type, Achene.twopt.angle.dev.max, pap.tip.cyp.max.area.ratio,
                               Pap.len.ach.len.ratio, pap.tip.ecc, coneSA, Plant_ID, Species)

dt_traits_merge_fre2$Species <- dplyr::recode(dt_traits_merge_fre2$Species, `L. fremontii` = "fre")

dt_traits_merge_gla2 <- dplyr::select(dt_traits_merge_gla, Velocity_Vals_m_sec, Seed_Mass_mg, 
                                 Seed_Type, Achene.twopt.angle.dev.max, pap.tip.cyp.max.area.ratio,
                                 Pap.len.ach.len.ratio, pap.tip.ecc, coneSA, Plant_ID, Species)

dt_traits_merge_gla2$Species <- dplyr::recode(dt_traits_merge_gla2$Species, `L. glaberrima` = "gla")


#' Make a composite dataframe with all species 
dt_traits_merge_cfg <- full_join(dt_traits_merge_cal2, dt_traits_merge_fre2) %>% full_join(dt_traits_merge_gla2) 

#' #### Seed mass (sm) plot 
#'

#' Calculate the marginal effect (slope) and the standard error in the slope for each species
sm_effect_cal <- as.data.frame(effect("Seed_Mass_mg", mod = cal_id_gamma_alt, xlevels=list(Seed_Mass_mg=seq(0.01,0.14,0.01))))
sm_effect_fre <- as.data.frame(effect("Seed_Mass_mg", mod = fre_id_gamma_alt, xlevels=list(Seed_Mass_mg=seq(0.00,0.10,0.01))))
sm_effect_gla <- as.data.frame(effect("Seed_Mass_mg", mod = gla_id_gamma_alt, xlevels=list(Seed_Mass_mg=seq(0.03,0.43,0.02))))
sm_effect_fre
0.7504345 -1.6798769
#' Add the species ID as a column
sm_effect_cal$Species <- "cal"
sm_effect_fre$Species <- "fre"
sm_effect_gla$Species <- "gla"

#' Create the composite data frame for plotting 
sm_effect <- full_join(sm_effect_cal, sm_effect_fre) %>% full_join(sm_effect_gla)

#' Calculate the percentage drop in flight distance per 1 unit difference (marginal effect size)
pcnt_eff_sm_cal <- mean(PcntEffectVec(sm_effect_cal, 10))
pcnt_eff_sm_cal
PcntEffectVec
pcnt_eff_sm_fre <- mean(PcntEffectVec(sm_effect_fre, 10))
pcnt_eff_sm_fre

pcnt_eff_sm_gla <- mean(PcntEffectVec(sm_effect_gla, 5))
pcnt_eff_sm_gla

#' Create the plot
sm_plot <- ggplot(data = sm_effect, aes(Seed_Mass_mg, fit, fill = Species)) +
  geom_ribbon(aes(ymin=fit-se, ymax=fit+se)) + #, alpha= 0.8
  geom_line(colour = "black") +
  geom_point(data = dt_traits_merge_cfg, aes(Seed_Mass_mg, Velocity_Vals_m_sec, color = Species, fill = Species)) + #, alpha = 0.4
  scale_colour_manual(values = c("#d95f02", "#21908C", "#440154")) +
  scale_fill_manual(values = c("#d95f02", "#21908C", "#440154")) + 
  scale_x_continuous(limits=c(0, 0.45), breaks = seq(0,0.45,0.05)) +
  scale_y_continuous(limits=c(0.45, 2.80), breaks = seq(0.5,2.80,0.25)) + 
  labs(
    #title = "2019 Lasthenia dispersal kernels",
    x = "Seed mass (mg)",
    y = "Drop velocity (m/s)"
  ) + 
  #geom_text(aes(x = 0.25, y = 1, label = "β = 8.60, p = 1.67e-10"), size = 3.5, colour = viridis(3)[3]) + #cal
  #geom_text(aes(x = 0.25, y = 0.75, label = "β = 9.80, p = 3.49e-06"), size = 3.5, colour = viridis(3)[2]) + #fre
  #geom_text(aes(x = 0.25, y = 0.5, label = "β = 2.62, p = 8.39e-04"), size = 3.5, colour = viridis(3)[1]) + #gla 
  guides(fill = F, colour = F) +
  theme_classic()

sm_plot

#' #### Angle of the Cypsela (ca) (deviation from straight = 180 degrees)
#' 
range(dt_traits_merge_cal$Achene.twopt.angle.dev.max, na.rm = TRUE)
#' Calculate the marginal effect (slope) and the standard error in the slope for each species
ca_effect_cal <- as.data.frame(effect("Achene.twopt.angle.dev.max", mod = cal_id_gamma_alt, xlevels=list(Achene.twopt.angle.dev.max=seq(0,26,0.5))))
ca_effect_fre <- as.data.frame(effect("Achene.twopt.angle.dev.max", mod = fre_id_gamma_alt, xlevels=list(Achene.twopt.angle.dev.max=seq(4.00,28.00,0.5))))
ca_effect_gla <- as.data.frame(effect("Achene.twopt.angle.dev.max", mod = gla_id_gamma_alt, xlevels=list(Achene.twopt.angle.dev.max=seq(1.00,16.10,0.5))))
ca_effect_gla
#' Add the species ID as a column
ca_effect_cal$Species <- "cal"
ca_effect_fre$Species <- "fre"
ca_effect_gla$Species <- "gla"

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
ca_plot <- ggplot(data = ca_effect, aes(Achene.twopt.angle.dev.max, fit, fill = Species)) +
  geom_ribbon(aes(ymin=fit-se, ymax=fit+se)) + #, alpha = 0.8
  geom_line(colour = "black") +
  geom_point(data = dt_traits_merge_cfg, aes(Achene.twopt.angle.dev.max, Velocity_Vals_m_sec, color = Species, fill = Species))+ #, alpha = 0.4
  scale_colour_manual(values = c("#d95f02", "#21908C", "#440154")) +
  scale_fill_manual(values = c("#d95f02", "#21908C", "#440154")) +
  scale_x_continuous(limits=c(0, 30), breaks = seq(0,30,5)) +
  scale_y_continuous(limits=c(0.45, 2.80), breaks = seq(0.5,2.80,0.25)) +  
  labs(
    #title = "2019 Lasthenia dispersal kernels",
    x = "Angle of the cypsela (°)",
    y = ""
  ) + 
  theme(axis.title.y = element_blank()) +
  #geom_text(aes(x = 10, y = 2.50, label = "β = -0.0100, p = 0.0310"), size = 3.5, hjust= 0, colour = viridis(3)[3]) + #cal
  #geom_text(aes(x = 10, y = 2.25, label = "β = -0.0105, p = 0.0647"), size = 3.5, hjust= 0, colour = viridis(3)[2]) + #fre
  #geom_text(aes(x = 10, y = 2, label = "β = 2.20e-03, p = 0.504"), size = 3.5, hjust= 0, colour = viridis(3)[1]) + #gla 
  theme_classic() +
  guides(fill = F, colour = F) 

ca_plot

#' #### Pappus area to cypsela area ratio
#' 
#' Calculate the marginal effect (slope) and the standard error in the slope for each species
ar_effect_cal <- as.data.frame(effect("pap.tip.cyp.max.area.ratio", mod = cal_id_gamma_alt, xlevels=list(pap.tip.cyp.max.area.ratio=seq(1,44,1))))
ar_effect_fre <- as.data.frame(effect("pap.tip.cyp.max.area.ratio", mod = fre_id_gamma_alt, xlevels=list(pap.tip.cyp.max.area.ratio=seq(0,13,0.5))))
ar_effect_gla <- as.data.frame(effect("pap.tip.cyp.max.area.ratio", mod = gla_id_gamma_alt, xlevels=list(pap.tip.cyp.max.area.ratio=seq(0,3,0.1))))
ar_effect_gla

summary(gla_id_gamma_alt)
confint(gla_id_gamma_alt, method = "Wald")

#' Add the species ID as a column
ar_effect_cal$Species <- "cal"
ar_effect_fre$Species <- "fre"
ar_effect_gla$Species <- "gla"

#' Calculate the percentage drop in flight distance per 1 unit difference (marginal effect size)
pcnt_eff_ar_cal <- mean(PcntEffectVec(ar_effect_cal, 1))
pcnt_eff_ar_cal

pcnt_eff_ar_fre <- mean(PcntEffectVec(ar_effect_fre, 2))
pcnt_eff_ar_fre

pcnt_eff_ar_gla <- mean(PcntEffectVec(ar_effect_gla, 10))
pcnt_eff_ar_gla

#' Create the composite data frame for plotting
ar_effect <- full_join(ar_effect_cal, ar_effect_fre) %>% full_join(ar_effect_gla)

#' Create the plot
ar_plot <- ggplot(data = ar_effect, aes(pap.tip.cyp.max.area.ratio, fit, fill = Species)) +
  geom_ribbon(aes(ymin=fit-se, ymax=fit+se)) + #, alpha = 0.8
  geom_line(colour="black")+
  geom_point(data = dt_traits_merge_cfg, aes(pap.tip.cyp.max.area.ratio, Velocity_Vals_m_sec, color = Species, fill = Species))+  #, alpha = 0.4
  scale_colour_manual(values = c("#d95f02", "#21908C", "#440154"), labels = c("L.californica","L. fremontii","L. glaberrima")) +
  scale_fill_manual(values = c("#d95f02", "#21908C", "#440154"), labels = c("L.californica","L. fremontii","L. glaberrima")) + 
  scale_x_continuous(limits=c(0, 45), breaks = seq(0,45,5)) +
  scale_y_continuous(limits=c(0.45, 2.80), breaks = seq(0.5,2.80,0.25)) + 
  labs(
    #title = "2019 Lasthenia dispersal kernels",
    x = "Pappus tip area (mm^2)/cypsela cs area (mm^2)",
    y = ""
  ) + 
  #geom_text(aes(x = 15, y = 2.25, label = "β = -3.54e-03, p = 0.321"), size = 3.5, hjust = 0, colour = viridis(3)[3]) + #cal
  #geom_text(aes(x = 15, y = 2.00, label = "β = 5.97e-03, p = 0.533"), size = 3.5, hjust = 0, colour = viridis(3)[2]) + #fre 
  #geom_text(aes(x = 15, y = 1.75, label = "β = 0.151, p = 0.0181"), size = 3.5, hjust = 0, colour = viridis(3)[1]) + #gla 
  theme_classic() + 
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(face = "italic"))

ar_plot

#' #### Pappus length to cypsela length ratio
#' 
#' Calculate the marginal effect (slope) and the standard error in the slope for each species
lr_effect_cal <- as.data.frame(effect("Pap.len.ach.len.ratio", mod = cal_id_gamma_alt, xlevels=list(Pap.len.ach.len.ratio=seq(0.8,1.8,0.1))))
lr_effect_fre <- as.data.frame(effect("Pap.len.ach.len.ratio", mod = fre_id_gamma_alt, xlevels=list(Pap.len.ach.len.ratio=seq(0.1,1.8,0.1))))
lr_effect_gla <- as.data.frame(effect("Pap.len.ach.len.ratio", mod = gla_id_gamma_alt, xlevels=list(Pap.len.ach.len.ratio=seq(0.3,0.7,0.1))))

#' Add the species ID as a column
lr_effect_cal$Species <- "cal"
lr_effect_fre$Species <- "fre"
lr_effect_gla$Species <- "gla"

#' Create the composite data frame for plotting
lr_effect <- full_join(lr_effect_cal, lr_effect_fre) %>% full_join(lr_effect_gla)

#' Create the plot
lr_plot <- ggplot(data = lr_effect, aes(Pap.len.ach.len.ratio, fit, fill = Species)) +
  geom_ribbon(aes(ymin=fit-se, ymax=fit+se)) + #, alpha = 0.8
  geom_line(colour = "black")+
  geom_point(data = dt_traits_merge_cfg, aes(Pap.len.ach.len.ratio, Velocity_Vals_m_sec, color = Species, fill = Species))+ #, alpha = 0.4
  scale_colour_manual(values = c("#d95f02", "#21908C", "#440154")) +
  scale_fill_manual(values = c("#d95f02", "#21908C", "#440154")) + 
  scale_y_continuous(limits=c(0.45, 2.80), breaks = seq(0.5,2.80,0.25)) + 
  scale_x_continuous(limits=c(0, 2), breaks = seq(0,2,0.5)) +
  labs(
    #title = "2019 Lasthenia dispersal kernels",
    x = "Pappus length (mm)/cypsela length (mm)",
    y = "Drop velocity (m/s)"
  ) + 
  #geom_text(aes(x = 0.75, y = 2.40, label = "β = -0.278, p = 0.131"), size = 3.5, colour = viridis(3)[3]) + #cal
  #geom_text(aes(x = 0.75, y = 2.25, label = "β = 0.0321, p = 0.808"), size = 3.5, colour = viridis(3)[2]) + #fre 
  #geom_text(aes(x = 0.75, y = 2.10, label = "β = 0.140, p = 0.849"), size = 3.5, colour = viridis(3)[1]) + #gla 
  guides(fill = F, colour = F)+
  theme_classic()

lr_plot

#' #### Variance in nflorescence height plot
#' 
#' Calculate the marginal effect (slope) and the standard error in the slope for each species
pe_effect_cal <- as.data.frame(effect("pap.tip.ecc", mod = cal_id_gamma_alt, xlevels=list(pap.tip.ecc=seq(0,1,0.01))))
pe_effect_fre <- as.data.frame(effect("pap.tip.ecc", mod = fre_id_gamma_alt, xlevels=list(pap.tip.ecc=seq(0,1,0.01))))
pe_effect_gla <- as.data.frame(effect("pap.tip.ecc", mod = gla_id_gamma_alt, xlevels=list(pap.tip.ecc=seq(0.1,1,0.01))))

#' Add the species ID as a column
pe_effect_cal$Species <- "cal"
pe_effect_fre$Species <- "fre"
pe_effect_gla$Species <- "gla"

#' Create the composite data frame for plotting
pe_effect <- full_join(pe_effect_cal, pe_effect_fre) %>% full_join(pe_effect_gla)

#' Create the plot
pe_plot <- ggplot(data = pe_effect, aes(pap.tip.ecc, fit, fill = Species)) +
  geom_ribbon(aes(ymin=fit-se, ymax=fit+se)) + #, alpha = 0.8
  geom_line(colour = "black")+
  geom_point(data = dt_traits_merge_cfg, aes(pap.tip.ecc, Velocity_Vals_m_sec, color = Species, fill = Species))+ #, alpha = 0.4
  scale_colour_manual(values = c("#d95f02", "#21908C", "#440154")) +
  scale_fill_manual(values = c("#d95f02", "#21908C", "#440154")) + 
  scale_y_continuous(limits=c(0.45, 2.80), breaks = seq(0.5,2.80,0.25)) +   
  scale_x_continuous(limits=c(0, 1.00), breaks = seq(0, 1, 0.25)) +
  labs(
    #title = "2019 Lasthenia dispersal kernels",
    x = "Pappus tip eccentricity",
    y = ""
  ) + 
  theme_classic() +
  #geom_text(aes(x = 0.50, y = 2.40, label = "β = 0.0226, p = 0.833"), size = 3.5, colour = viridis(3)[3]) + 
  #geom_text(aes(x = 0.50, y = 2.25, label = "β = 0.0380, p = 0.724"), size = 3.5, colour = viridis(3)[2]) + 
  #geom_text(aes(x = 0.50, y = 2.10, label = "β = -0.0999, p = 0.536"), size = 3.5, colour = viridis(3)[1]) + 
  guides(fill = F, colour = F)

pe_plot

#' #### Surface area of the cypsela
#' 
#' Calculate the marginal effect (slope) and the standard error in the slope for each species
sa_effect_cal <- as.data.frame(effect("coneSA", mod = cal_id_gamma_alt, xlevels=list(coneSA=seq(0.3,0.9,0.01))))
sa_effect_fre <- as.data.frame(effect("coneSA", mod = fre_id_gamma_alt, xlevels=list(coneSA=seq(0.2,1.0,0.01))))
sa_effect_gla <- as.data.frame(effect("coneSA", mod = gla_id_gamma_alt, xlevels=list(coneSA=seq(0.3,2.06,0.01))))

#' Add the species ID as a column
sa_effect_cal$Species <- "cal"
sa_effect_fre$Species <- "fre"
sa_effect_gla$Species <- "gla"

#' Create the composite data frame for plotting
sa_effect <- full_join(sa_effect_cal, sa_effect_fre) %>% full_join(sa_effect_gla)

#' Create the plot
sa_plot <- ggplot(data = sa_effect, aes(coneSA, fit, fill = Species)) +
  geom_ribbon(aes(ymin=fit-se, ymax=fit+se)) + #, alpha = 0.8
  geom_line(colour = "black")+
  geom_point(data = dt_traits_merge_cfg, aes(coneSA, Velocity_Vals_m_sec, color = Species, fill = Species))+ #, alpha = 0.4
  scale_colour_manual(values = c("#d95f02", "#21908C", "#440154")) +
  scale_fill_manual(values = c("#d95f02", "#21908C", "#440154")) +  
  scale_y_continuous(limits=c(0.45, 2.80), breaks = seq(0.5,2.80,0.25)) + 
  scale_x_continuous(limits=c(0.25,2.25), breaks = seq(0.25,2.25,0.25)) +
  labs(
    #title = "2019 Lasthenia dispersal kernels",
    x = "Cypsela surface area (mm^2)",
    y = ""
  ) + 
  theme_classic() +
  #geom_text(aes(x = 1.25, y = 2.40, label = "β = -0.133, p = 0.701"), size = 3.5, colour = viridis(3)[3]) + 
  #geom_text(aes(x = 1.25, y = 2.25, label = "β = 1.52e-03, p = 0.997"), size = 3.5, colour = viridis(3)[2]) + 
  #geom_text(aes(x = 1.25, y = 2.10, label = "β = 0.0498, p = 0.771"), size = 3.5, colour = viridis(3)[1]) + 
  guides(fill = F, colour = F)

sa_plot

summary(gla_id_gamma_alt)

#' Use the patchwrok package to make this plot

layout <- c(
  area(1, 2, 1, 5),
  area(2, 1, 3, 2),
  area(2, 3, 3, 4),
  area(2, 5, 3, 6),
  area(4, 1, 5, 2),
  area(4, 3, 5, 4),
  area(4, 5, 5, 6))

composite_plot <- guide_area() + lr_plot + ar_plot + pe_plot + sm_plot + ca_plot + sa_plot + plot_layout(design = layout, guides = "collect")

composite_plot

ggsave(plot = composite_plot, "/Users/Courtney/Documents/Thesis Chapter 1/Spatial dispersal paper/Figures/DT_scatter_composite_plot.eps", width = 15, height = 11)


vignette("models", "emmeans")
