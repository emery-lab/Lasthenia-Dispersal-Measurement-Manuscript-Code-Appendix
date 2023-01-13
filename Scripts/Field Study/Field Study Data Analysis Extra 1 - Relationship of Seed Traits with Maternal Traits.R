#### Field Study Data Analysis Extra 1: Relationship of Seed Traits with Maternal Traits

# Courtney Van Den Elzen

# Initiated: July 2020
# Latest Update: January 2023 (syntax update for Jan 2023 package version compatibility)

### ------------------
### Load Packages -----

library(tidyverse)
library(lme4)
library(car)
library(MASS)

### ------------------
### Read In Data ----

# Seed weight data for the neighboring plants
seed_mass_data <- read_csv("./Data/Field Study/Field Study Cleaned Data - Seed Mass vs Maternal Plant Size Relationship Data.csv")

# Flower count data for the neighboring plants
flwr_count_data <- 
     read_csv("./Data/2019 Field Data Jepson Neighbour Seed Counts.csv") %>%
     dplyr::select(-c("...17")) %>%
     dplyr::filter(!is.na(flwr_count)) %>%
     dplyr::group_by(plant_id) %>%
     dplyr::summarize(tot_flwrs = sum(flwr_count))


# Plant length data for the neighboring plants
plant_length_data <- read_csv("./Data/Field Study Cleaned Data - Neighbouring Plant Data.csv")

# L. californica
modsel_cal <- read_csv("./Data/Field Study/modsel_cal_analysis_pt5.csv")

# L. glaberrima
modsel_gla <- read_csv("./Data/Field Study/modsel_gla_analysis_pt5.csv")



