#### Field Dispersal Data Analysis Extra 4: Predicting Proportions of Seed Recovered 3

# Courtney Van Den Elzen

# Initiated: October 2019
# Latest Update: January 2023 (syntax update for Jan 2023 package version compatibility)

# Description: Use chosen predictive model to predict the number of viable seeds 
#              in each 2019 focal maternal plant

# WARNING: You must first run Field Study Data Analysis Extra 3 - Predicting 
#          Proportions of Seed Recovered.R so the code has access to the models


### ------------------
### Load Packages -----
library(ciTools)
library(foreach)
library(doParallel)
library(tidyverse)
library(lme4)

### ------------------
### Custom functions ----- 

## Bootstrap the prediction interval for the poisson glm
boot_pi <- function(model, pdata, n, p) {
  odata <- model$data
  lp <- (1 - p) / 2
  up <- 1 - lp
  set.seed(2016)
  seeds <- round(runif(n, 1, 1000), 0)
  boot_y <- foreach(i = 1:n, .combine = rbind) %dopar% {
    set.seed(seeds[i])
    bdata <- odata[sample(seq(nrow(odata)), size = nrow(odata), replace = TRUE), ]
    bpred <- predict(update(model, data = bdata), type = "response", newdata = pdata)
    rpois(length(bpred), lambda = bpred)
  }
  boot_ci <- t(apply(boot_y, 2, quantile, c(lp, up)))
  return(data.frame(pred = predict(model, newdata = pdata, type = "response"), lower = boot_ci[, 1], upper = boot_ci[, 2]))
}


### ------------------
### Read In Combined Data (from Part 1) -----

# Read in the diameter data for the maternal plants
mat_plants <- read_csv("./Data/Field Study/Field Study Raw Data - Maternal Inflorescences.csv")

# Tidy the data
mat_plants_good <- mat_plants %>%
  dplyr::mutate(Species = case_when(Species == "L. californica" ~ "cal", 
                             Species == "L. fremontii" ~ "fre",
                             Species == "L. glaberrima" ~ "gla"),
         Size_num = as.numeric(Size)) %>%
  dplyr::mutate(Inflor_ID = paste(Hub, Species,"T", Transect, "M", Mom, "I", Inflorescence, sep = ""),
         Mat_ID = paste(Hub, Species,"T", Transect, "M", Mom, sep = ""))


### ------------------
### Predict Seed Numbers -------

# Create a dataframe of the predictor variables
mat_diams_spp <- data.frame(Species = mat_plants_good$Species, 
                            Diameter = mat_plants_good$Size,
                            Mat_ID = mat_plants_good$Mat_ID) %>%
  as_tibble() %>%
  dplyr::mutate(Diameter = as.numeric(Diameter)) %>%
  dplyr::filter(!is.na(Diameter))


## ------
## Poisson GLM -------
pois_pred_w_PI <- boot_pi(disc_poisson, mat_diams_spp, 1000, 0.95)


## ------
## Negative Binomial GLM -------

mat_disc_pred_nb <- predict(disc_nb, newdata = mat_diams_spp, type = "response")

nb_pred_df <- data.frame(Disc_tube_count = mat_disc_pred_nb, 
                         Species = mat_diams_spp$Species, 
                         Diameter = mat_diams_spp$Diameter)

nb_pred_w_PI <- ciTools::add_pi(nb_pred_df, 
                                disc_nb, 
                                names = c("lpb", "upb")) %>% 
  dplyr::mutate(Disc_tube_count = NULL, Species = NULL, Diameter = NULL)


## ------
## Negative Gaussian LM -------

gauss_pred_w_PI <- predict(disc_gauss, newdata = mat_diams_spp, type = "response",  interval = "prediction")


## ------
## Finalizing Predicted Seed Numbers -------

# Moving forward with the poisson model since it had lower prediction error than 
# the nb model. The gaussian model had the lowest MSE but it also has prediction 
# intervals that include negative values so I'm going to stick with the poisson.

mat_inflor_data <- 
  mat_diams_spp %>% 
  dplyr::mutate(pois_nseed_pred = pois_pred_w_PI$pred, 
                pois_nseed_pred_lower = pois_pred_w_PI$lower, 
                pois_nseed_pred_upper = pois_pred_w_PI$upper)


