#### 2019 Field Dispersal Data Analysis 5: Statistical Modeling of Combined Data

# Courtney Van Den Elzen

# Initiated: December 2019
# Latest Update: January 2023 (syntax update for Jan 2023 version compatibility)

# Description: Statistical model fitting and comparison for all field study data


### ------------------
### Load Packages -------
library(tidyverse)
library(lme4)
library(nlme)
library(car)
library(effects)


### ------------------
### Read In Combined Data (all cleaning and joining done in scripts 1-4) -------

kern_hw_sv_mp_midf <- read_csv("./Data/Field Study/Field Study Cleaned Data - Combined Kernel, Maternal, Sur Veg, and Wind Speed Data.csv")

### ------------------
### Clean Data for Modeling -------

kern_hw_sv_mp_midf <- kern_hw_sv_mp_midf %>% 
  mutate(
    SppID = factor(SppID), 
    Full_ID = factor(Full_ID), 
    Hub = factor(Hub),
    Transect = factor(Transect),
    Mom = factor(Mom),
    Pin = factor(Pin)
  )

# Change data types for use in modeling functions
kern_hw_sv_mp_midf$SppID <- as.factor(kern_hw_sv_mp_midf$SppID)

# Change the reference level (needed to do all pairwise statistical comparisons by species)
kern_hw_sv_mp_midf_frebase <- within(kern_hw_sv_mp_midf, SppID <- relevel(SppID, ref = "fre"))

# Break the data down by species
modsel_cal <- filter(kern_hw_sv_mp_midf, SppID == "cal")
modsel_fre <- filter(kern_hw_sv_mp_midf, SppID == "fre")
modsel_gla <- filter(kern_hw_sv_mp_midf, SppID == "gla")

### ------------------
### Data Modeling: Dispersal Distance By Species -------

# Description: fit models of dispersal distance (coded as Distance) or inter-seed
#              spread (coded as ptcentdist) as the response variable and species 
#              as the predictor. Models are Gamma-distributed GLMMs with identity, 
#              log, or inverse link functions. The same model is fit with 
#              L. californica or L. fremontii as the reference level for the 
#              species variable in order to get p-values for all pairwise species 
#              comparisons


## ------
## Distance - identity link, L. californica base -------

# Gamma GLMM model (identity link) of species dispersal distance differences, 
# controlling for plant id
spp_mod_id <- glmer(data = kern_hw_sv_mp_midf, 
                    Distance ~ SppID + (1|Full_ID), 
                    family = Gamma(link = "identity"))

# summary of model fit
summary(spp_mod_id)

# Wald confidence intervals for all fitted parameter estimates
confint(spp_mod_id, method = "Wald")

# overall significance of the species predictor
anova(spp_mod_id)

# marginal predictor effects from this model using the 'effects' package
spp_mod_id_eff <- effects::effect("SppID", spp_mod_id, se = T)
spp_mod_id_eff # fitted values
spp_mod_id_eff$se # standard errors


## ------
## Distance - identity link, L. fremontii base -------

# Gamma GLMM model (identity link) of species dispersal distance differences, 
# controlling for plant id, L. fremontii as reference level
spp_mod_id_fb <- glmer(data = kern_hw_sv_mp_midf_frebase, 
                       Distance ~ SppID + (1|Full_ID), 
                       family = Gamma(link = "identity"))

# summary of model fit
summary(spp_mod_id_fb)

# Wald confidence intervals for all fitted parameter estimates
confint(spp_mod_id_fb, method = "Wald")

# effects from this model. 
spp_mod_id_fb_eff <- effects::effect("SppID", spp_mod_id_fb, se = TRUE)
spp_mod_id_fb_eff # fitted values
spp_mod_id_fb_eff$se # standard errors


## ------
## Inter-seed spread - identity link, L. californica base -------

# Model of inter-seed species by species, controlling for plant id
spp_mod_id_disp <- glmer(data = kern_hw_sv_mp_midf, 
                         ptcentdist ~ SppID + (1|Full_ID), 
                         family = Gamma(link = "identity"))

# update function used to extend model fitting procedure (helps with convergence)
spp_mod_id_disp <- update(spp_mod_id_disp, 
                          start=getME(spp_mod_id_disp, c("theta","fixef")), 
                          control=glmerControl(optimizer="bobyqa", 
                                               optCtrl=list(maxfun=5e6)))

# model fit summary
summary(spp_mod_id_disp)

# Wald confidence intervals for all fitted parameter estimates
confint(spp_mod_id_disp, method = "Wald")

# effects from this model. 
spp_mod_id_disp_eff <- effects::effect("SppID", spp_mod_id_disp, se = TRUE)
spp_mod_id_disp_eff # fitted values
spp_mod_id_disp_eff$se # standard errors


## ------
## Inter-seed spread - identity link, L. fremontii base -------

# Again but with L. fremontii as the baseline (to get p-values for each contrast)
spp_mod_id_disp_fb <- glmer(data = kern_hw_sv_mp_midf_frebase, 
                            ptcentdist ~ SppID + (1|Full_ID), 
                            family = Gamma(link = "identity"))

# convergence (see above)
spp_mod_id_disp_fb <- update(spp_mod_id_disp_fb, 
                             start=getME(spp_mod_id_disp_fb, c("theta","fixef")), 
                             control=glmerControl(optimizer="bobyqa", 
                                                  optCtrl=list(maxfun=5e6)))
# model fit summary
summary(spp_mod_id_disp_fb)

# Wald confidence intervals for all fitted parameter estimates
confint(spp_mod_id_disp_fb, method = "Wald")

#' marginal effects from this model. 
spp_mod_id_disp_fb_eff <- effects::effect("SppID", spp_mod_id_disp, se = TRUE)
spp_mod_id_disp_fb_eff # fitted values
spp_mod_id_disp_fb_eff$se # standard errors

### ------------------
### Data Modeling: L. californica Relationships of Traits with Dispersal Distance -------

# Description: fit models of dispersal distance (coded as Distance) as the 
#              response variable and maternal traits and environmental 
#              covariates as the predictors. Models are Gamma-distributed GLMMs 
#              with identity, log, or inverse link functions.

## ------
## Distance - identity link -------
modsel_cal$mean_surveg_length_cm
modsel_cal$mean_prop_speed_pcnt

# Identity link
cal_id <- glmer(Distance ~ mean_prop_speed_pcnt + mean_surveg_length_cm + mean_inflor_height_cm + n_inflor + (1|Full_ID), 
                    family = Gamma(link = "identity"), 
                    control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=6e8)),
                    data = modsel_cal)

# model summary
summary(cal_id)

# update helps with model fit convergence
cal_id <- update(cal_id, start=getME(cal_id, c("theta","fixef")), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=8e8)))

## ------
## Distance - log link -------

cal_log <- glmer(Distance ~ mean_prop_speed_pcnt + mean_surveg_length_cm + mean_inflor_height_cm + n_inflor + (1|Full_ID), 
                 family = Gamma(link = "log"),
                 control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=6e8)), 
                 data = modsel_cal)

cal_log <- update(cal_log, start=getME(cal_log, c("theta","fixef")), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=8e8)))

## ------
## Distance - inverse link -------
cal_inv <- glmer(Distance ~ mean_prop_speed_pcnt + mean_surveg_length_cm + mean_inflor_height_cm + n_inflor + (1|Full_ID), 
                 family = Gamma(link = "inverse"), 
                 control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=6e8)),
                 data = modsel_cal)

cal_inv <- update(cal_inv, start=getME(cal_inv, c("theta","fixef")), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=6e8)))


## -----
## Model comparisons using AIC -------

#' Log link and ID link models are best. 
AIC(cal_id)
AIC(cal_log)
AIC(cal_inv)

#' Summary of the gamma glmm with identity link for *L. californica*
summary(cal_id)
confint(cal_id, method = "Wald")


### ------------------
### Data Modeling: L.  fremontii Relationships of Traits with Dispersal Distance -------

# Description: fit models of dispersal distance (coded as Distance) or inter-seed
#              spread (coded as ptcentdist) as the response variable and maternal
#              traits and environmental covariates as the predictors. Models are
#              Gamma-distributed GLMMs with identity, log, or inverse link 
#              functions.

## ------
## Distance - identity link -------

fre_id <- glmer(Distance ~ mean_prop_speed_pcnt + mean_surveg_length_cm + mean_inflor_height_cm + (1|Full_ID), 
                family = Gamma(link = "identity"),
                control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=6e8)),
                data = modsel_fre)


## ------
## Distance - log link -------

fre_log <- glmer(Distance ~  mean_prop_speed_pcnt + mean_surveg_length_cm + mean_inflor_height_cm + (1|Full_ID), 
                 family = Gamma(link = "log"), 
                 control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=6e8)),
                 data = modsel_fre)

fre_log <- update(fre_log, start=getME(fre_log, c("theta","fixef")), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=6e8)))


## ------
## Distance - inverse link -------

fre_inv <- glmer(Distance ~ mean_prop_speed_pcnt + mean_surveg_length_cm + mean_inflor_height_cm + (1|Full_ID), 
                 family = Gamma(link = "inverse"), 
                 control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=6e8)),
                 data = modsel_fre)

fre_inv <- update(fre_inv, 
                  start=getME(fre_inv, c("theta","fixef")), 
                  control=glmerControl(optimizer="bobyqa", 
                                       optCtrl=list(maxfun=6e8)))


## -----
## Model comparisons using AIC -------

AIC(fre_id) 
AIC(fre_log) 
AIC(fre_inv)

# Summary of the gamma glmm with identity link for L. fremontii
summary(fre_id)

confint(fre_id, method = "Wald")


### ------------------
### Data Modeling: L.  glaberrima Relationships of Traits with Dispersal Distance -------

# Description: fit models of dispersal distance (coded as Distance) or inter-seed
#              spread (coded as ptcentdist) as the response variable and maternal
#              traits and environmental covariates as the predictors. Models are
#              Gamma-distributed GLMMs with identity, log, or inverse link 
#              functions.


## ------
## Distance - identity link -------

gla_id <- glmer(Distance ~ mean_prop_speed_pcnt + mean_surveg_length_cm + mean_inflor_height_cm + var_inflor_height_cm + n_inflor + (1|Full_ID), 
                    family = Gamma(link = "identity"),
                    data = modsel_gla)

gla_id <- update(gla_id, 
                 start=getME(gla_id, c("theta","fixef")), 
                 control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=6e8)))


## ------
## Distance - log link -------

gla_log <- glmer(Distance ~ mean_prop_speed_pcnt + mean_surveg_length_cm + n_inflor + mean_inflor_height_cm + var_inflor_height_cm + (1|Full_ID), 
                     family = Gamma(link = "log"), 
                     control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=6e8)),
                     data = modsel_gla)

gla_log <- update(gla_log, start=getME(gla_log, c("theta","fixef")), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=6e8)))


## ------
## Distance - inverse link -------
gla_inv <- glmer(Distance ~  mean_prop_speed_pcnt + mean_surveg_length_cm + mean_inflor_height_cm + var_inflor_height_cm + n_inflor + (1|Full_ID), 
                 family = Gamma(link = "inverse"), 
                 data = modsel_gla)

gla_inv <- update(gla_inv, start=getME(gla_inv, c("theta","fixef")), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=6e8)))


## ------
## Model comparisons using AIC -------

# ID link is best
AIC(gla_id)
AIC(gla_log)
AIC(gla_inv)

# Summary of the gamma glmm with log link for *L. glaberrima*
summary(gla_id)

confint(gla_id, method = "Wald")

# Rerun the model without the outlier for variance
modsel_gla_filt <- filter(modsel_gla, mean_surveg_length_cm < 35)

gla_id_filt <- glmer(Distance ~ mean_prop_speed_pcnt + mean_surveg_length_cm + mean_inflor_height_cm + var_inflor_height_cm + n_inflor + (1|Full_ID), 
                family = Gamma(link = "identity"),
                data = modsel_gla_filt)

# model summary
summary(gla_id_filt)

### ------------------
### Data Modeling: L. californica relationships of Traits with Inter-Seed Distance -------

# Description: fit models of inter-seed spread (coded as ptcentdist) as the 
#              response variable and maternal traits and environmental 
#              covariates as the predictors. Models are Gamma-distributed GLMMs 
#              with identity, log, or inverse link functions.

## ------
## Inter-seed Spread - identity link ------
cal_id_disp <- glmer(ptcentdist ~ mean_prop_speed_pcnt + mean_surveg_length_cm + mean_inflor_height_cm + n_inflor + (1|Full_ID), 
                     family = Gamma(link = "identity"), 
                     control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=8e8)),
                     data = modsel_cal)

## ------
## Inter-seed Spread - log link ------
cal_log_disp <- glmer(ptcentdist ~ mean_prop_speed_pcnt + mean_surveg_length_cm + mean_inflor_height_cm + n_inflor + (1|Full_ID), 
                          family = Gamma(link = "log"),
                          control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=6e8)),
                          data = modsel_cal)

cal_log_disp <- update(cal_log_disp, start=getME(cal_log_disp, c("theta","fixef")), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=6e8)))


## ------
## Inter-seed Spread - inverse link ------
cal_inv_disp <- glmer(ptcentdist ~ mean_prop_speed_pcnt + mean_surveg_length_cm + mean_inflor_height_cm + n_inflor + (1|Full_ID), 
                      family = Gamma(link = "inverse"), 
                      control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=8e8)),
                      data = modsel_cal)

## ------
## Model comparisons using AIC ------
AIC(cal_id_disp)
AIC(cal_log_disp) 
AIC(cal_inv_disp)

# Summary of the gamma glmm with id link for L. californica
summary(cal_id_disp)

### ------------------
### Data Modeling: L. fremontii relationships of Traits with Inter-Seed Distance -------

# Description: fit models of inter-seed spread (coded as ptcentdist) as the 
#              response variable and maternal traits and environmental 
#              covariates as the predictors. Models are Gamma-distributed GLMMs 
#              with identity, log, or inverse link functions.


## ------
## Inter-seed Spread - identity link ------
fre_id_disp <- glmer(ptcentdist ~  mean_prop_speed_pcnt + mean_surveg_length_cm + mean_inflor_height_cm + (1|Full_ID), 
                     family = Gamma(link = "identity"),
                     control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=6e8)),
                     data = modsel_fre)


## ------
## Inter-seed Spread - log link ------
fre_log_disp <- glmer(ptcentdist ~  mean_prop_speed_pcnt + mean_surveg_length_cm + mean_inflor_height_cm + n_inflor + (1|Full_ID), 
                      family = Gamma(link = "log"),
                      control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=6e8)),
                      data = modsel_fre)


## ------
## Inter-seed Spread - inverse link ------
fre_inv_disp <- glmer(ptcentdist ~  mean_prop_speed_pcnt + mean_surveg_length_cm + mean_inflor_height_cm + n_inflor +(1|Full_ID), 
                      family = Gamma(link = "inverse"),
                      control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=6e8)),
                      data = modsel_fre)


## ------
## Model comparisons using AIC ------
AIC(fre_id_disp) 
AIC(fre_log_disp)
AIC(fre_inv_disp) 

# Summary of the gamma glmm with identity link for L. fremontii
summary(fre_id_disp)


### ------------------
### Data Modeling: L. glaberrima relationships of Traits with Inter-Seed Distance -------

# Description: fit models of inter-seed spread (coded as ptcentdist) as the 
#              response variable and maternal traits and environmental 
#              covariates as the predictors. Models are Gamma-distributed GLMMs 
#              with identity, log, or inverse link functions.


## ------
## Inter-seed Spread - identity link ------
gla_id_disp <- glmer(ptcentdist ~  mean_prop_speed_pcnt + mean_surveg_length_cm + mean_inflor_height_cm + var_inflor_height_cm + n_inflor + (1|Full_ID), 
                     family = Gamma(link = "identity"), 
                     #control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=6e8)),
                     data = modsel_gla)

## ------
## Inter-seed Spread - log link ------
gla_log_disp <- glmer(ptcentdist ~ mean_prop_speed_pcnt + mean_surveg_length_cm + mean_inflor_height_cm + var_inflor_height_cm + n_inflor + (1|Full_ID), 
                      family = Gamma(link = "log"),
                      control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=6e8)),
                      data = modsel_gla)

## ------
## Inter-seed Spread - inverse link ------
gla_inv_disp <- glmer(ptcentdist ~  mean_prop_speed_pcnt + mean_surveg_length_cm + mean_inflor_height_cm + var_inflor_height_cm + n_inflor + (1|Full_ID), 
                      family = Gamma(link = "inverse"), 
                      control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=6e8)),
                      data = modsel_gla)

## ------
## Model comparisons using AIC ------
AIC(gla_id_disp)
AIC(gla_log_disp)
AIC(gla_inv_disp)

# Summary of the gamma glmm with log link for *L. glaberrima*
summary(gla_id_disp)

# Rerun the model without the outlier for variance
modsel_gla_filt <- filter(modsel_gla, var_inflor_height_cm < 9)

gla_id_filt_disp <- glmer(ptcentdist ~ mean_prop_speed_pcnt + mean_surveg_length_cm + mean_inflor_height_cm + var_inflor_height_cm + (1|Full_ID), 
                     family = Gamma(link = "identity"), 
                     control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=6e8)),
                     data = modsel_gla_filt)

gla_id_filt_disp <- update(gla_id_filt_disp, start=getME(gla_id_filt_disp, c("theta","fixef")), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=7e8)))


summary(gla_id_filt_disp)

anova(gla_id_filt_disp,gla_id_filt_disp)

## -------------
## Write Data to File 

#write_csv(modsel_cal, "./Data/Field Study/Field Study Cleaned Data - Californica Species Modeling Data Subset.csv")
#write_csv(modsel_cal, "./Data/Field Study/Field Study Cleaned Data - Fremontii Species Modeling Data Subset.csv")
#write_csv(modsel_cal, "./Data/Field Study/Field Study Cleaned Data - Glaberrima Species Modeling Data Subset.csv")
