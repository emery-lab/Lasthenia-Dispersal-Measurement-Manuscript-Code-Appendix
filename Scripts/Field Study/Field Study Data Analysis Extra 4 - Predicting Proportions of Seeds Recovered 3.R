#' # Field Study 2019 Jepson Prairie
#' # Recovered Seed Proportions Part 3
#' 
#' ### Part 3/4: Use chosen predictive model to predict the number of viable seeds in each 2019 focal maternal plant
#' 
#' Courtney Van Den Elzen
#' 
#' October 21st 2019 - Nov 22nd 2019
#' 
#' Part 2 must be run first
#' 
#' Libraries
library(foreach)
library(doParallel)

#' #### Custom functions

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

#' #### Read in and tidy the data
#'
#' Read in the diameter data for the maternal plants
mat_plants <- read_csv("/Users/Courtney/Documents/Thesis Chapter 1/Field Study/2019 Full Study Files/2019 Field Data Jepson Maternal Inflorescence Traits.csv")

#' Tidy the data
mat_plants_good <- mat_plants %>% filter(Omit_from_analysis == "no") %>% 
  mutate(Species = case_when(Species == "L. californica" ~ "cal", 
                             Species == "L. fremontii" ~ "fre",
                             Species == "L. glaberrima" ~ "gla"),
         Size = as.numeric(Size)) %>%
  mutate(Inflor_ID = paste(Hub, Species,"T", Transect, "M", Mom, "I", Inflorescence, sep = ""),
         Mat_ID = paste(Hub, Species,"T", Transect, "M", Mom, sep = ""),
         Year = as.character("2019"))

#' #### Predict maternal plant seed numbers using each of the fitted models
#' 
#' Create a dataframe of the predictor variables
mat_diams_spp <- data.frame(Species = mat_plants_good$Species, Diameter = mat_plants_good$Size, Year = mat_plants_good$Year)

head(mat_diams_spp)

#' poisson model (load boot_pi function below)
pois_pred_w_PI <- boot_pi(disc_poisson, mat_diams_spp, 1000, 0.95)

head(pois_pred_w_PI)

#' negative binomial model
mat_disc_pred_nb <- predict(disc_nb, newdata = mat_diams_spp, type = "response")
nb_pred_df <- data.frame(Disc_tube_count = mat_disc_pred_nb, Species = mat_diams_spp$Species, Diameter = mat_diams_spp$Diameter, Year = mat_diams_spp$Year)
nb_pred_w_PI <- ciTools::add_pi(nb_pred_df, disc_nb, names = c("lpb", "upb")) %>% mutate(Disc_tube_count = NULL, Species = NULL, Diameter = NULL, Year = NULL)

head(nb_pred_w_PI)

#' gaussian model
gauss_pred_w_PI <- predict(disc_gauss, newdata = mat_diams_spp, type = "response",  interval = "prediction")

head(gauss_pred_w_PI)

#' Moving forward with the poisson model since it had lower prediction error than the nb model. 
#' The gaussian model had the lowest MSE but it also has prediction intervals that include negative values so 
#' I'm going to stick with the poisson.
mat_inflor_data <- mat_plants_good %>% 
  mutate(pois_nseed_pred = round(pois_pred_w_PI$pred), 
         pois_nseed_pred_lower = pois_pred_w_PI$lower, 
         pois_nseed_pred_upper = pois_pred_w_PI$upper) %>%
  dplyr::select(-Notes)



#' For knitting to html: knitr::spin("/Users/Courtney/Documents/Thesis Chapter 1/Field Study/2019 Full Study Files/Scripts/Predicting Found Proportions 3.R")
