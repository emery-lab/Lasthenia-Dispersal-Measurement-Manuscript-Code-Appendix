#### 2019 Field Dispersal Data Analysis Extra 3: Predicting Proportions of Seed Recovered 2

# Courtney Van Den Elzen

# Initiated: October 2019
# Latest Update: January 2023 (syntax update for Jan 2023 package version compatibility)

# Description: Model the disc tube counts and choose best predictive model using
#              drop-one cross validation


### ------------------
### Load Packages -----
library(lmtest)

### ------------------
### Read In Combined Data (from Part 1) -----

nseed_data_all_good <- 
  read_csv("./Data/Field Study/Field Study Cleaned Data - 2018 and 2019 Neighbour Data Cleaned and Combined.csv") 

nseed_data_all_good$Species

### ------------------
### Split Data by Species -----

# Split by species for modeling
nseed_data_all_good_cal <- dplyr::filter(nseed_data_all_good, Species == "cal")
nseed_data_all_good_fre <- dplyr::filter(nseed_data_all_good, Species == "fre")
nseed_data_all_good_gla <- dplyr::filter(nseed_data_all_good, Species == "gla")

### ------------------
### Data Modeling -------

# Model the predictions of viable disc seed number (will not be modeling ray seeds b/c they don't disperse)

## ------
## Poisson GLM -------

disc_poisson <- glm(Disc_tube_count ~ Species*Diameter*Year,
                        data = nseed_data_all_good,
                        family = "poisson")

## ------
## Negative Binomial GLM -------
disc_nb <- glm.nb(Disc_tube_count ~ Species*Diameter*Year, data = nseed_data_all_good)


## ------
## Regular (Gaussian) LM -------
disc_gauss <- lm(Disc_tube_count ~ Species*Diameter*Year, data = nseed_data_all_good)


### ------------------
### Prediction Error Rate -------

# Measuring prediction error rate between the two models - drop 1 cross validation

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


# Poisson glm model MSE for disc seeds
pois_err

# Negative binomial glm model MSE for disc seeds
nb_err

# Gaussian linear model MSE for disc seeds
gauss_err

