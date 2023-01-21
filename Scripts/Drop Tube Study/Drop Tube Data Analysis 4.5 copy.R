#' # Drop Tube Data Analysis Part 4
#' 
#' ### Visualization of seed trait effects and modeling trait effects on mean velocity
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


summarize(group_by(dt_traits_merge, Plant_ID), n())
dt_traits_merge

#' #### Load data (run part 3 first) and clean
dt_traits_merge_cal <- dt_traits_merge %>% filter(Species == "L. californica")
dt_traits_merge_fre <- dt_traits_merge %>% filter(Species == "L. fremontii")
dt_traits_merge_gla <- dt_traits_merge %>% filter(Species == "L. glaberrima")

#' Create a column that is a unique ID of the parent plant to use as a random effect
dt_traits_merge_cal <- dt_traits_merge_cal %>% mutate(Plant_ID = glue("{Set}cal{Number}m{Mom}n{Neighbour}"))
dt_traits_merge_fre <- dt_traits_merge_fre %>% mutate(Plant_ID = glue("{Set}fre{Number}m{Mom}n{Neighbour}"))
dt_traits_merge_gla <- dt_traits_merge_gla %>% mutate(Plant_ID = glue("{Set}gla{Number}m{Mom}n{Neighbour}"))

#' #### Load data (run part 3 first) and clean
dt_traits_merge_cal_disc <- dt_traits_merge %>% filter(Species == "L. californica") %>% filter(Seed_Type == "d")
dt_traits_merge_fre_disc <- dt_traits_merge %>% filter(Species == "L. fremontii") %>% filter(Seed_Type == "d")
dt_traits_merge_gla_disc <- dt_traits_merge %>% filter(Species == "L. glaberrima") %>% filter(Seed_Type == "d")

#' Create a column that is a unique ID of the parent plant to use as a random effect
dt_traits_merge_cal_disc <- dt_traits_merge_cal_disc %>% mutate(Plant_ID = glue("{Set}cal{Number}m{Mom}n{Neighbour}"))
dt_traits_merge_fre_disc <- dt_traits_merge_fre_disc %>% mutate(Plant_ID = glue("{Set}fre{Number}m{Mom}n{Neighbour}"))
dt_traits_merge_gla_disc <- dt_traits_merge_gla_disc %>% mutate(Plant_ID = glue("{Set}gla{Number}m{Mom}n{Neighbour}"))

#' Merge all species datasets
dt_traits_merge_2 <- bind_rows(dt_traits_merge_cal_disc, dt_traits_merge_fre_disc, dt_traits_merge_gla_disc)
dt_traits_merge_2


dim(dt_traits_merge_cal_disc)

#write.csv(dt_traits_merge_2, "/Users/Courtney/Documents/Thesis Chapter 1/Drop Tube Study/Data/Drop Tube Velocity Traits All Species.csv")

#' ### Modeling
#' 
#' Run a gamma GLMM on californica (Identity link function had the lowest AIC score)
#' 
#' Model with identity link
cal_id_gamma <- glmer(Velocity_Vals_m_sec ~ Seed_Mass_mg + Seed_Type + Achene.twopt.angle.dev.max + pap.tip.cyp.max.area.ratio + 
                        Pap.len.ach.len.ratio + pap.tip.ecc + coneSA + (1|Plant_ID), 
                      data = dt_traits_merge_cal, 
                      family=Gamma(link="identity"), 
                      control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=5e8)))

#' Model with log link
cal_log_gamma <- glmer(Velocity_Vals_m_sec ~ Seed_Mass_mg + Achene.twopt.angle.dev.max + pap.tip.cyp.max.area.ratio + 
                         Pap.len.ach.len.ratio + pap.tip.ecc + coneSA + Seed_Type + (1|Plant_ID), 
                       data = dt_traits_merge_cal, 
                       family=Gamma(link="log"), 
                       control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=5e8)))

#' Model with inverse link
cal_inv_gamma <- glmer(Velocity_Vals_m_sec ~ Seed_Mass_mg + Achene.twopt.angle.dev.max + pap.tip.cyp.max.area.ratio + 
                         Pap.len.ach.len.ratio + pap.tip.ecc + coneSA + Seed_Type + (1|Plant_ID), 
                       data = dt_traits_merge_cal, 
                       family=Gamma(link="inverse"), 
                       control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=5e8)))

#' Model with identity link
cal_id_gamma_alt <- glmer(Velocity_Vals_m_sec ~ Seed_Mass_mg + Achene.twopt.angle.dev.max + pap.tip.cyp.max.area.ratio + 
                          Pap.len.ach.len.ratio + pap.tip.ecc + coneSA + (1|Plant_ID), 
                        data = dt_traits_merge_cal_disc, 
                        family=Gamma(link="identity"),
                        control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=5e8)))

#' Model with log link
cal_log_gamma_alt <- glmer(Velocity_Vals_m_sec ~ Seed_Mass_mg + Achene.twopt.angle.dev.max + pap.tip.cyp.max.area.ratio + 
                           Pap.len.ach.len.ratio + pap.tip.ecc + coneSA + (1|Plant_ID), 
                         data = dt_traits_merge_cal, 
                         family=Gamma(link="log"), 
                         control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=5e8)))

#' Model with inverse link
cal_inv_gamma_alt <- glmer(Velocity_Vals_m_sec ~ Seed_Mass_mg + Achene.twopt.angle.dev.max + pap.tip.cyp.max.area.ratio + 
                           Pap.len.ach.len.ratio + pap.tip.ecc + coneSA + (1|Plant_ID), 
                         data = dt_traits_merge_cal, 
                         family=Gamma(link="inverse"), 
                         control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=5e8)))

AIC(cal_id_gamma) # ***
AIC(cal_log_gamma)
AIC(cal_inv_gamma)

AIC(cal_id_gamma_alt) 
AIC(cal_log_gamma_alt)
AIC(cal_inv_gamma_alt)

summary(cal_id_gamma)
summary(cal_id_gamma_alt)

#' This is weird, but the shape parameter is 1/sigma^2 (see this link here from Ben Bolker: https://stat.ethz.ch/pipermail/r-sig-mixed-models/2015q1/023375.html)
#' 
#' all of these seem to match up with the data pretty well (plotted on Wolfram Alpha) 
shape_cal_id_gamma <- 1/sigma(cal_id_gamma)^2
disp_cal_id_gamma <- 1/shape_cal_id_gamma
scale_cal_id_gamma <- sigma(cal_id_gamma) / sqrt(shape_cal_id_gamma)

#' Run a gamma GLM on fremontii (log link had the lowest AIC score (inverse wouldn't converge))
#' 
#' Model with identity link
fre_id_gamma <- glmer(Velocity_Vals_m_sec ~ Seed_Mass_mg + Achene.twopt.angle.dev.max + pap.tip.cyp.max.area.ratio + 
                        Pap.len.ach.len.ratio + pap.tip.ecc + coneSA + Seed_Type + (1|Plant_ID), 
                      data = dt_traits_merge_fre, 
                      family=Gamma(link="identity"), 
                      control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=6e8)))

#' Model with long link
fre_log_gamma <- glmer(Velocity_Vals_m_sec ~ Seed_Mass_mg + Achene.twopt.angle.dev.max + pap.tip.cyp.max.area.ratio + 
                         Pap.len.ach.len.ratio + pap.tip.ecc + coneSA + Seed_Type + (1|Plant_ID), 
                       data = dt_traits_merge_fre, 
                       family=Gamma(link="log"), 
                       control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=6e8)))

#' Model with identity link
fre_id_gamma_alt <- glmer(Velocity_Vals_m_sec ~ Seed_Mass_mg + Achene.twopt.angle.dev.max + pap.tip.cyp.max.area.ratio + 
                        Pap.len.ach.len.ratio + pap.tip.ecc + coneSA + (1|Plant_ID), 
                      data = dt_traits_merge_fre_disc, 
                      family=Gamma(link="identity"), 
                      control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=6e8)))

range(dt_traits_merge_fre_disc$Seed_Mass_mg)
range(dt_traits_merge_fre_disc$Achene.twopt.angle.dev.max)
range(dt_traits_merge_fre_disc$pap.tip.cyp.max.area.ratio)
range(dt_traits_merge_fre_disc$Pap.len.ach.len.ratio)
range(dt_traits_merge_fre_disc$pap.tip.ecc)
range(dt_traits_merge_fre_disc$coneSA)

fre_id_gamma_alt <- update(fre_id_gamma_alt, start=getME(fre_id_gamma_alt, c("theta","fixef")), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=6e8)))


#' Model with long link
fre_log_gamma_alt <- glmer(Velocity_Vals_m_sec ~ Seed_Mass_mg + Achene.twopt.angle.dev.max + pap.tip.cyp.max.area.ratio + 
                         Pap.len.ach.len.ratio + pap.tip.ecc + coneSA + (1|Plant_ID), 
                       data = dt_traits_merge_fre_disc, 
                       family=Gamma(link="log"), 
                       control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=6e8)))


AIC(fre_id_gamma)
AIC(fre_log_gamma)

AIC(fre_id_gamma_alt)
AIC(fre_log_gamma_alt)

summary(fre_id_gamma)
summary(fre_id_gamma_alt)

#' Run a gamma GLM on glaberrima
#' 
#' Chose log link as it had the lowest AIC score (identity link would not run)
gla_id_gamma <- glmer(Velocity_Vals_m_sec ~ Seed_Mass_mg + Achene.twopt.angle.dev.max + 
                        pap.tip.cyp.max.area.ratio + Pap.len.ach.len.ratio + pap.tip.ecc + 
                        coneSA + (1|Plant_ID), 
                      data = dt_traits_merge_gla,
                      family=Gamma(link="identity"),
                      control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=6e8)))

gla_id_gamma <- update(gla_id_gamma, start=getME(gla_id_gamma, c("theta","fixef")), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=6e8)))

gla_log_gamma <- glmer(Velocity_Vals_m_sec ~ Seed_Mass_mg + Achene.twopt.angle.dev.max + pap.tip.cyp.max.area.ratio +
                         Pap.len.ach.len.ratio + pap.tip.ecc + coneSA + (1|Plant_ID), 
                       data = dt_traits_merge_gla, 
                       family=Gamma(link="log"),
                       control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=6e8)))

gla_inv_gamma <- glmer(Velocity_Vals_m_sec ~ Seed_Mass_mg + Achene.twopt.angle.dev.max + pap.tip.cyp.max.area.ratio +
                         Pap.len.ach.len.ratio + pap.tip.ecc + coneSA + (1|Plant_ID), 
                       data = dt_traits_merge_gla, 
                       family=Gamma(link="inverse"),
                       control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=6e8)))

AIC(gla_id_gamma)
AIC(gla_log_gamma)
AIC(gla_inv_gamma)


#' #### Showing the "winning" species models
#' 
#' L. californica
summary(cal_id_gamma)
#' L. fremontii
summary(fre_id_gamma)
#' L. glaberrima
summary(gla_id_gamma)

#' #### Checking for multicollinearity
#' 
#' I was concerned about possibly multicollinearity in the predictors for this analysis, so here I constructed ordinary least squares regression models for each of the species and then calculated the variance inflation factors for each of the predictors. 
#' 
#' Raw correlations between the variables
#'
#' L. cal - nothing of concern
cor(drop_na(dt_traits_merge_cal[c("Seed_Mass_mg", "Achene.twopt.angle.dev.max", "pap.tip.cyp.max.area.ratio", "Pap.len.ach.len.ratio", "pap.tip.ecc", "coneSA")]))

#' L. fre - nothing of concern
cor(drop_na(dt_traits_merge_fre[c("Seed_Mass_mg", "Achene.twopt.angle.dev.max", "pap.tip.cyp.max.area.ratio", "Pap.len.ach.len.ratio", "pap.tip.ecc", "coneSA")]))

#' L. gla - nothing of concern except 0.70 cor between coneSA and seed mass. 
cor(drop_na(dt_traits_merge_gla[c("Seed_Mass_mg", "Achene.twopt.angle.dev.max", "pap.tip.cyp.max.area.ratio", "Pap.len.ach.len.ratio", "pap.tip.ecc", "coneSA")]))

#' Calculating VIF between parameters
#'
#' L. cal OLS regression
cal_lm <- lm(Velocity_Vals_m_sec ~ Seed_Mass_mg + Achene.twopt.angle.dev.max + pap.tip.cyp.max.area.ratio +
               Pap.len.ach.len.ratio + pap.tip.ecc, 
             data = dt_traits_merge_cal)

#' VIF for the cal model - nothing of concern
vif(cal_lm)

#' L. fre OLS regression
fre_lm <- lm(Velocity_Vals_m_sec ~ Seed_Mass_mg + Achene.twopt.angle.dev.max + pap.tip.cyp.max.area.ratio +
               Pap.len.ach.len.ratio + pap.tip.ecc, 
             data = dt_traits_merge_fre)

#' VIF for the fre model - nothing of concern
vif(fre_lm)

#' L. gla OLS regression
gla_lm <- lm(Velocity_Vals_m_sec ~ Seed_Mass_mg + Achene.twopt.angle.dev.max + pap.tip.cyp.max.area.ratio +
               Pap.len.ach.len.ratio + pap.tip.ecc, 
             data = dt_traits_merge_gla)

#' VIF for the gla model - nothing of concern
vif(gla_lm)

#knitr::spin("/Users/Courtney/Documents/Thesis Chapter 1/Drop Tube Study/Analysis/Drop Tube Data Analysis 4.R")
