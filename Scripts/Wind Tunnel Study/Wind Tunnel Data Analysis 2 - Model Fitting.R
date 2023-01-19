#### Wind Tunnel Data Analysis 2: Trait Data Cleaning and Aggregation

# Courtney Van Den Elzen

# Original: 2018 

# Latest Update: January 2023 (syntax update for Jan 2023 version compatibility)

# Description: This script explores the trait data through plots, builds Gamma 
#              GLMM models of distance traveled in the wind tunnel using seed 
#              traits as predictors

# WARNING: You must run Wind Tunnel Data Analysis 1 first to have access to the data

# The seed traits used are:
 
# 1. Seed mass in mg ("Seed.weight")
# 2. Achene curvature (i.e. angle) measured as max deviance from a straight line ("Achene.twopt.angle.dev.max")
# 3. Ratio of the area of the tip of pappus to the max cross-sectional area of the achene ("pap.tip.cyp.max.area.ratio")
# 4. Ratio of the pappus length to the achene length ("Pap.len.ach.len.ratio")
# 5. Eccentricity (ellipsoidalness) of the pappus tip ("pap.tip.ecc")
# 6. Surface area of the achene ("coneSA")


### ------------------
### Load Packages -------

library(cubature)
library(lme4)
library(tidyverse)
library(viridis)


### ------------------
### Split and Summarize Data -------

# Look at the data (cleaned in previous script)
colnames(traits_mod_agg.cal.mod) # Trait data
colnames(dists2_main.cal) # Distance and identity data

# Subset the data to the traits I care about
traits_mod_agg.cal.mod2 <- traits_mod_agg.cal.mod[,c("Tube.ID", "Species", "Wild.Grnhs", "Seed.weight", "Achene.twopt.angle.dev.max", "pap.tip.cyp.max.area.ratio", "pap.tip.pap.base.area.ratio", "Pap.perp.max.len.mean", "Ach.len.lin.mean", "Pap.len.ach.len.ratio", "pap.tip.ecc", "coneSA")] 
traits_mod_agg.fre.mod2 <- traits_mod_agg.fre.mod[,c("Tube.ID", "Species", "Wild.Grnhs", "Seed.weight", "Achene.twopt.angle.dev.max", "pap.tip.cyp.max.area.ratio", "pap.tip.pap.base.area.ratio", "Pap.perp.max.len.mean", "Ach.len.lin.mean", "Pap.len.ach.len.ratio", "pap.tip.ecc", "coneSA")] 
traits_mod_agg.gla.mod2 <- traits_mod_agg.gla.mod[,c("Tube.ID", "Species", "Wild.Grnhs", "Seed.weight", "Achene.twopt.angle.dev.max", "pap.tip.cyp.max.area.ratio", "pap.tip.pap.base.area.ratio", "Pap.perp.max.len.mean", "Ach.len.lin.mean", "Pap.len.ach.len.ratio", "pap.tip.ecc", "coneSA")] 

# Merge the dataframes
dist_traits_cal <- merge(traits_mod_agg.cal.mod2, dists2_main.cal, by = "Tube.ID")
dist_traits_fre <- merge(traits_mod_agg.fre.mod2, dists2_main.fre, by = "Tube.ID")
dist_traits_gla <- merge(traits_mod_agg.gla.mod2, dists2_main.gla, by = "Tube.ID")

# Filter to only greenhouse-grown seeds
dist_traits_cal_grnhs <- dplyr::filter(dist_traits_cal, Wild.Grnhs == "Greenhouse")
dist_traits_fre_grnhs <- dplyr::filter(dist_traits_fre, Wild.Grnhs == "Greenhouse")
dist_traits_gla_grnhs <- dplyr::filter(dist_traits_gla, Wild.Grnhs == "Greenhouse")

# Filter to only wild-grown seeds
dist_traits_cal_wild <- dplyr::filter(dist_traits_cal, Wild.Grnhs == "Wild")
dist_traits_fre_wild <- dplyr::filter(dist_traits_fre, Wild.Grnhs == "Wild")
dist_traits_gla_wild <- dplyr::filter(dist_traits_gla, Wild.Grnhs == "Wild")

# Join by species
dist_traits_all <- full_join(dist_traits_cal, dist_traits_fre) %>% full_join(dist_traits_gla)


### ------------------
### Data Plots -------

## -----------------
## L. californica Plots ------

# Achene.twopt.angle.dev.max by plant for L. californica
ggplot(data = dist_traits_cal_wild, aes(x = Plant.ID, y = Achene.twopt.angle.dev.max, color = Plant.ID)) + 
  geom_point() + 
  geom_line() + 
  ggtitle("L. californica achene angle deviation wild seeds")

# Pap.len.ach.len.ratio by plant ID for L. californica
ggplot(data = dist_traits_cal_wild, aes(x = Plant.ID, y = Pap.len.ach.len.ratio, color = Plant.ID)) + 
  geom_point() +
  geom_line() +
  ggtitle("L. californica pappus length/achene length wild seeds")

# Seed.weight by plant ID for L. californica
ggplot(data = dist_traits_cal_wild, aes(x = Plant.ID, y = Seed.weight, color = Plant.ID)) + 
  geom_point()+
  geom_line() + 
  ggtitle("L. californica seed mass wild seeds")

# Achene.twopt.angle.dev.max by plant ID for L. californica
ggplot(data = dist_traits_cal_wild, aes(x = Plant.ID, y = Achene.twopt.angle.dev.max, color = Plant.ID)) + 
  geom_point()+
  geom_line() + 
  ggtitle("L. californica achene curvature")

# pap.tip.ecc by plant ID for L. californica
ggplot(data = dist_traits_cal_wild, aes(x = Plant.ID, y = pap.tip.ecc, color = Plant.ID)) + 
  geom_point()+
  geom_line() + 
  ggtitle("L. californica pappus eccentricity")

# coneSA by plant ID for L. californica
ggplot(data = dist_traits_cal_wild, aes(x = Plant.ID, y = coneSA, color = Plant.ID)) + 
  geom_point()+
  geom_line() + 
  ggtitle("L. californica achene surface area")


## -----------------
## L. fremontii Plots ------

# Achene.twopt.angle.dev.max by plant for L. fremontii
ggplot(data = dist_traits_fre_wild, aes(x = Plant.ID, y = Achene.twopt.angle.dev.max, color = Plant.ID)) + 
  geom_point() + 
  geom_line() + 
  ggtitle("L. fremontii achene angle deviation wild seeds")

# Pap.len.ach.len.ratio by plant ID for L. fremontii
ggplot(data = dist_traits_fre_wild, aes(x = Plant.ID, y = Pap.len.ach.len.ratio, color = Plant.ID)) + 
  geom_point()+
  geom_line() + 
  ggtitle("L. fremontii pappus length/achene length wild seeds")

# Seed.weight by plant ID for L. californica
ggplot(data = dist_traits_fre_wild, aes(x = Plant.ID, y = Seed.weight, color = Plant.ID)) + 
  geom_point()+
  geom_line() +
  ggtitle("L. fremontii seed mass wild seeds")

# Achene.twopt.angle.dev.max by plant ID for L. fremontii
ggplot(data = dist_traits_fre_wild, aes(x = Plant.ID, y = Achene.twopt.angle.dev.max, color = Plant.ID)) + 
  geom_point()+
  geom_line() + 
  ggtitle("L. fremontii achene curvature")

# pap.tip.ecc by plant ID for L. fremontii
ggplot(data = dist_traits_fre_wild, aes(x = Plant.ID, y = pap.tip.ecc, color = Plant.ID)) + 
  geom_point()+
  geom_line() + 
  ggtitle("L. fremontii pappus eccentricity")

# coneSA by plant ID for L. fremontii
ggplot(data = dist_traits_fre_wild, aes(x = Plant.ID, y = coneSA, color = Plant.ID)) + 
  geom_point()+
  geom_line() + 
  ggtitle("L. fremontii achene surface area")


## -----------------
## L. glaberrima Plots ------

# Achene.twopt.angle.dev.max by plant for L. glaberrima
ggplot(data = dist_traits_gla_wild, aes(x = Plant.ID, y = Achene.twopt.angle.dev.max, color = Plant.ID)) + 
  geom_point() + 
  geom_line() + 
  ggtitle("L. glaberrima achene angle deviation wild seeds")

# Pap.len.ach.len.ratio by plant ID for L. glaberrima
ggplot(data = dist_traits_gla_wild, aes(x = Plant.ID, y = Pap.len.ach.len.ratio, color = Plant.ID)) + 
  geom_point()+
  geom_line() + 
  ggtitle("L. glaberrima pappus length/achene length wild seeds")

# Seed.weight by plant ID for L. californica
ggplot(data = dist_traits_gla_wild, aes(x = Plant.ID, y = Seed.weight, color = Plant.ID)) + 
  geom_point()+
  geom_line() +
  ggtitle("L. glaberrima seed mass wild seeds")

# Achene.twopt.angle.dev.max by plant ID for L. glaberrima
ggplot(data = dist_traits_gla_wild, aes(x = Plant.ID, y = Achene.twopt.angle.dev.max, color = Plant.ID)) + 
  geom_point()+
  geom_line() + 
  ggtitle("L. glaberrima achene curvature")

# pap.tip.ecc by plant ID for L. glaberrima
ggplot(data = dist_traits_gla_wild, aes(x = Plant.ID, y = pap.tip.ecc, color = Plant.ID)) + 
  geom_point()+
  geom_line() + 
  ggtitle("L. glaberrima pappus eccentricity")

# coneSA by plant ID for L. glaberrima
ggplot(data = dist_traits_gla_wild, aes(x = Plant.ID, y = coneSA, color = Plant.ID)) + 
  geom_point() +
  geom_line() + 
  ggtitle("L. glaberrima achene surface area")

### ------------------
### Model Fitting -------

## -----------------
## L. californica Models ------

# inverse link
calmod_inv <- glm(Total.Hypoteneuse.Distance ~ Seed.weight + Achene.twopt.angle.dev.max +
                   pap.tip.cyp.max.area.ratio + Pap.len.ach.len.ratio +
                   pap.tip.ecc + coneSA,
                 data = dist_traits_cal, 
                 family=Gamma(link="inverse"))

# identity link
calmod_id <- glm(Total.Hypoteneuse.Distance ~ Seed.weight + Achene.twopt.angle.dev.max +
                         pap.tip.cyp.max.area.ratio + Pap.len.ach.len.ratio +
                         pap.tip.ecc + coneSA,
                       data = dist_traits_cal, 
                       family=Gamma(link="identity"))

# log link
calmod_log <- glm(Total.Hypoteneuse.Distance ~ Seed.weight + Achene.twopt.angle.dev.max +
                          pap.tip.cyp.max.area.ratio + Pap.len.ach.len.ratio +
                          pap.tip.ecc + coneSA,
                        data = dist_traits_cal, 
                        family=Gamma(link="log"))



# In all models, log link is slightly favored
AIC(calmod_inv)
AIC(calmod_id)
AIC(calmod_log)

#' Summaries of each model
summary(calmod_id)
summary(calmod_log)


## -----------------
## L. fremontii Models ------

# inverse link
fremod_inv <- glm(Total.Hypoteneuse.Distance ~ Seed.weight + Achene.twopt.angle.dev.max +
                   pap.tip.cyp.max.area.ratio + Pap.len.ach.len.ratio +
                   pap.tip.ecc + coneSA,
                 data = dist_traits_fre,
                 family=Gamma(link="inverse"))


# identity link
fremod_id <- glm(Total.Hypoteneuse.Distance ~ Seed.weight + Achene.twopt.angle.dev.max +
                            pap.tip.cyp.max.area.ratio + Pap.len.ach.len.ratio +
                            pap.tip.ecc + coneSA,
                          data = dist_traits_fre, 
                          family=Gamma(link="identity"))

# log link
fremod_log <- glm(Total.Hypoteneuse.Distance ~ Seed.weight + Achene.twopt.angle.dev.max +
                             pap.tip.cyp.max.area.ratio + Pap.len.ach.len.ratio +
                             pap.tip.ecc + coneSA,
                           data = dist_traits_fre, 
                           family=Gamma(link="log"))

# In all models, log link is slightly favored
AIC(fremod_inv)
AIC(fremod_id)
AIC(fremod_log)

# Summaries of each model
summary(fremod_log)
summary(fremod_id)


## -----------------
## L. glaberrima Models ------

# inverse link
glamod_inv <- glm(Total.Hypoteneuse.Distance ~ Seed.weight + Achene.twopt.angle.dev.max +
                   pap.tip.cyp.max.area.ratio + Pap.len.ach.len.ratio +
                   pap.tip.ecc + coneSA + Wild.Grnhs,
                 data = dist_traits_gla, 
                 family=Gamma(link="inverse"))

# identity link
glamod_id <- glm(Total.Hypoteneuse.Distance ~ Seed.weight + Achene.twopt.angle.dev.max +
                            pap.tip.cyp.max.area.ratio + Pap.len.ach.len.ratio +
                            pap.tip.ecc + coneSA + Wild.Grnhs,
                          data = dist_traits_gla, 
                          family=Gamma(link="identity"))

# log link
glamod_log <- glm(Total.Hypoteneuse.Distance ~ Seed.weight + Achene.twopt.angle.dev.max +
                             pap.tip.cyp.max.area.ratio + Pap.len.ach.len.ratio +
                             pap.tip.ecc + coneSA + Wild.Grnhs,
                           data = dist_traits_gla, 
                           family=Gamma(link="log"))

# In all models, log link is slightly favored
AIC(glamod_inv)
AIC(glamod_id)
AIC(glamod_log)

# Summary of the model
summary(glamod_log)
summary(glamod_id)
