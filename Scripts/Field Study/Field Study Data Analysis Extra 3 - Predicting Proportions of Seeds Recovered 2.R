#' # Field Study 2019 Jepson Prairie
#' # Recovered Seed Proportions Part 2
#' 
#' ### Part 2/4: Model the disc tube counts and choose best predictive model using drop-one cross validation
#' 
#' Courtney Van Den Elzen
#' 
#' October 21st 2019 - Nov 22nd 2019
#' 
#' Libraries 

library(lmtest)
#' Part 1 must be run first.

#' #### Split by species for modeling
nseed_data_all_good_cal <- dplyr::filter(nseed_data_all_good, Species == "cal")
nseed_data_all_good_fre <- dplyr::filter(nseed_data_all_good, Species == "fre")
nseed_data_all_good_gla <- dplyr::filter(nseed_data_all_good, Species == "gla")

#' #### Modeling first pass
#' Model the predictions of viable disc seed number (will not be modeling ray seeds b/c they don't disperse)
#' 
#' Poisson glm
disc_poisson <- glm(Disc_tube_count ~ Species*Diameter*Year,
                        data = nseed_data_all_good,
                        family = "poisson")

# *** previously did a likelihood ratio test to check whether adding year into the model increased predictive accuracy and it did.

#' Negative binomial glm
disc_nb <- glm.nb(Disc_tube_count ~ Species*Diameter*Year, data = nseed_data_all_good)

# *** also did a LR test here. Same result.

#' Gaussian lm
disc_gauss <- lm(Disc_tube_count ~ Species*Diameter*Year, data = nseed_data_all_good)
summary(disc_gauss)

# *** also did a LR test here. Same result.

#' #### Prediction Error Rate 
#' Measuring prediction error rate between the two models - drop 1 cross validation

# Number of observations total
n_obs <- dim(nseed_data_all_good)[1]

# Number of observations to drop
ndrop <- 1

# Initialize matrices to store predicted data
drop_pred.pois <- matrix(nrow = n_obs, ncol = ndrop)
drop_pred.nb <- matrix(nrow = n_obs, ncol = ndrop)
drop_pred.gauss <- matrix(nrow = n_obs, ncol = ndrop)

# Initialize matrix to store observed data
drop_obs <- matrix(nrow = n_obs, ncol = ndrop)

# Iterate through all data points dropping ndrop at a time, then fit model and predict dropped pts
for (i in 1:n_obs){
  if (ndrop == 1){ #drop one data point at a time
    rand_drop <- i
    
    #training dataset
    train_set <- nseed_data_all_good[-rand_drop,]
    
    #testing dataset
    drop_set <- nseed_data_all_good[rand_drop,]
    
    #fit both model types (poisson and negative binomial)
    disc_poisson <- glm(Disc_tube_count ~ Species*Diameter*Year, data = train_set, family = "poisson")
    disc_nb <- glm.nb(Disc_tube_count ~ Species*Diameter*Year, data = train_set)
    disc_gauss <- lm(Disc_tube_count ~ Species*Diameter*Year, data = train_set)

    #predict the dropped values and store them
    dropdf <- data.frame(Species = drop_set$Species, Diameter = drop_set$Diameter, Year = drop_set$Year)
    drop_pred.pois[i,] <- predict(disc_poisson, dropdf, type = "response")
    drop_pred.nb[i,] <- predict(disc_nb, dropdf, type = "response")
    drop_pred.gauss[i,] <- predict(disc_gauss, dropdf, type = "response")
    drop_obs[i,] <- drop_set$Disc_tube_count
    
  } else { #more than one point dropped at a time
    rand_drop <- round(runif(ndrop, min=1, max=n_obs))
    
    #all remaining data points are training data set
    train_set <- nseed_data_all_good[-rand_drop,]
    drop_set <- nseed_data_all_good[rand_drop,]
    
    #fit both model types (poisson and negative binomial)
    disc_poisson <- glm(Disc_tube_count ~ Species*Diameter*Year, data = train_set, family = "poisson")
    disc_nb <- glm.nb(Disc_tube_count ~ Species*Diameter*Year, data = train_set)
    disc_gauss <- lm(Disc_tube_count ~ Species*Diameter*Year, data = train_set)
    
    #predict the dropped values and store them
    dropdf <- data.frame(Species = drop_set$Species, Diameter = drop_set$Diameter, Year = drop_set$Year)
    drop_pred.pois[i,] <- predict(disc_poisson, dropdf, type = "response")
    drop_pred.nb[i,] <- predict(disc_nb, dropdf, type = "response")
    drop_pred.gauss[i,] <- predict(disc_gauss, dropdf, type = "response")

    drop_obs[i,] <- drop_set$Disc_tube_count
  }
}

# calculate the MSE (CV score) for each model
pois_err <- sum((drop_obs - drop_pred.pois)^2)/n_obs
nb_err <- sum((drop_obs - drop_pred.nb)^2)/n_obs
gauss_err <- sum((drop_obs - drop_pred.gauss)^2)/n_obs

#' Lowest score means lowest prediction error. Gaussian model works best here. 
#Poisson glm model MSE for disc seeds
pois_err

#Negative binomial glm model MSE for disc seeds
nb_err

#Gaussian linear model MSE for disc seeds
gauss_err


#' For knitting to html: knitr::spin("/Users/Courtney/Documents/Thesis Chapter 1/Field Study/2019 Full Study Files/Scripts/Predicting Found Proportions 2.R")

