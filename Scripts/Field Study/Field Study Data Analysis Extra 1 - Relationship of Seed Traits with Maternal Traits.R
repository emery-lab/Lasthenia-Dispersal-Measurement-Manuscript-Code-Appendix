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
     read_csv("./Data/Field Study/Field Study Raw Data - Neighbour Seed Counts.csv") %>%
     dplyr::select(-c("...17")) %>%
     dplyr::filter(!is.na(flwr_count)) %>%
     dplyr::group_by(plant_id) %>%
     dplyr::summarize(tot_flwrs = sum(flwr_count))


# Plant length data for the neighboring plants
plant_length_data <- read_csv("./Data/Field Study/Field Study Cleaned Data - Neighbouring Plant Data.csv")

# L. californica
modsel_cal <- read_csv("./Data/Field Study/Field Study Cleaned Data - Californica Species Modeling Data Subset.csv")

# L. glaberrima
modsel_gla <- read_csv("./Data/Field Study/Field Study Cleaned Data - Glaberrima Species Modeling Data Subset.csv")

### ------------------
### Clean Neighbor Data -----

# Correct the neighbor data types to fit the analyses

# Neighbour seed mass data
seed_mass_data2 <- seed_mass_data %>% 
        mutate_at(c("orig_species_label", "true_species", "hub", "transect", "mom", 
                    "neighbour", "PlantID", "count type", "chosen for germination trial"), 
                  as.factor) %>%
        mutate_at(c("total seed mass"), as.numeric)

# Neighbour total branch length data
plant_length_data2 <- plant_length_data %>% 
        mutate_at(c("species", "hub", "transect", "mom", "neighbour", "PlantID"), as.factor) %>%
        mutate_at(c("Total branch length"), as.numeric) %>%
        rename(total_branch_length = `Total branch length`)


# Get rid of NAs in the data

seed_mass_data3 <- dplyr::filter(seed_mass_data2, !is.na(`total seed mass`))

plant_length_data3 <- dplyr::filter(plant_length_data2, !is.na(total_branch_length))


# Summarize the data for each plant in N1 (combine first flwr and other flwr data), N2 data already combined

# group by plantID and sum rows of seed mass and seed number for each plant
seed_mass_data_byplant <- group_by(seed_mass_data3, PlantID, true_species) %>% 
        dplyr::summarize(total_seed_mass = sum(`total seed mass`), 
                         total_seed_count = sum(`seed count`))

# get rid of unwanted ID columns
plant_length_data4 <- dplyr::select(plant_length_data3, -c("species", "hub", "transect", "mom", "neighbour"))

# join the cleaned seed mass and plant length dataframes, then filter out bad rows, and drop unused levels of true_species factor
joined_data_filtered <- dplyr::left_join(seed_mass_data_byplant, plant_length_data4) %>%
        dplyr::filter(true_species %in% c("L. californica", "L. fremontii", "L. glaberrima")) %>%
        droplevels()

# round to 3 decimal places
joined_data_filtered$mean_seed_mass <- round(joined_data_filtered$total_seed_mass/joined_data_filtered$total_seed_count,3)

# join the flower count data
joined_data_filtered2 <- 
        joined_data_filtered %>%
        dplyr::left_join(flwr_count_data, by = c("PlantID" = "plant_id"))

# Focal plant data

# Add a column to the focal plant data for the relative infloresence height variable for L. californica (used below)

# Mean inflor height / mean surrounding veg height
modsel_cal$rel_inflor_height_to_surveg_height <- modsel_cal$mean_inflor_height_cm/modsel_cal$mean_surveg_length_cm
modsel_gla$rel_inflor_height_to_surveg_height <- modsel_gla$mean_inflor_height_cm/modsel_gla$mean_surveg_length_cm


# Filter to only one line per plant for plant-level trait analysis (picked pin 1 arbitrarily, does not affect the analysis)

# L. californica
modsel_cal_plantlvl <- filter(modsel_cal, Pin == 1) %>% 
        filter(!is.na(n_inflor)) %>% 
        filter(!is.na(tot_mat_length_cm))

# L. glaberrima
modsel_gla_plantlvl <- filter(modsel_gla, Pin == 1) %>%
        filter(!is.na(n_inflor)) %>% 
        filter(!is.na(tot_mat_length_cm))

# Separate the joined data by species

joined_data_filtered_cal <- filter(joined_data_filtered2, true_species == "L. californica")
joined_data_filtered_gla <- filter(joined_data_filtered2, true_species == "L. glaberrima")

### ------------------
### Modeling: L. californica -----

## -------------------
## Do larger plants produce more inflorescences? -----

# Scatterplots -----

ggplot(modsel_cal_plantlvl, aes(x = tot_mat_length_cm, y = n_inflor)) + 
        geom_point()

# Modeling -----

# model
inflor.totlen.cal.lm <- lm(n_inflor ~ tot_mat_length_cm, data = modsel_cal_plantlvl)

# summary
summary(inflor.totlen.cal.lm)

# residual plots - wonky
plot(inflor.totlen.cal.lm)


## -------------------
## Do plants with more inflorescences grow in taller surrounding vegetation? -----
# Scatterplots -----

ggplot(modsel_cal_plantlvl, aes(x = mean_surveg_length_cm, y = n_inflor)) + 
        geom_point() 


# Modeling -----

# model
inflor.surveglen.cal.lm <- lm(n_inflor ~ mean_surveg_length_cm, data = modsel_cal_plantlvl)

# summary
summary(inflor.surveglen.cal.lm)

# residual plots - wonky
plot(inflor.surveglen.cal.lm)


## -------------------
# Do larger plants produce heavier seeds? -----

# Scatterplots -----

ggplot(joined_data_filtered_cal, aes(y = total_branch_length, x = mean_seed_mass)) + 
        geom_point()

# Modeling -----

seedmass.totlen.cal.lm <- lm(mean_seed_mass ~ total_branch_length, data = joined_data_filtered_cal)
summary(seedmass.totlen.cal.lm)
plot(seedmass.totlen.cal.lm )

# -------------------
# Do taller plants grow in taller vegetation? -----


# Scatterplots -----

# raw heights of both focal plants and surrounding vegetation length 
ggplot(modsel_cal_plantlvl, aes(x = mean_inflor_height_cm, y = mean_surveg_length_cm)) + 
        geom_point() + 
        geom_smooth(method= "lm")

# Relative height of focal plant relative to inflorenscence height 
ggplot(modsel_cal_plantlvl, aes(y = rel_inflor_height_to_surveg_height, x = mean_inflor_height_cm)) + 
        geom_point() + 
        geom_smooth(method= "lm")

# Modeling -----

# taller plants tend to occur in taller vegetation (absolute measurements)
mih.msvl.cal.lm <- lm(mean_inflor_height_cm ~ mean_surveg_length_cm, data = modsel_cal_plantlvl)

summary(mih.msvl.cal.lm)

# shorter plants (absolute height) tend to grow in relatively taller surrounding vegetation
relihsvh.msvl.cal.lm <- lm(rel_inflor_height_to_surveg_height ~ mean_inflor_height_cm, data = modsel_cal_plantlvl)

summary(relihsvh.msvl.cal.lm)

# -------------------
# Do larger plants grow in taller vegetation? -----

# Scatterplots -----

ggplot(modsel_cal_plantlvl, aes(x = tot_mat_length_cm, y = mean_surveg_length_cm)) +
        geom_point() + 
        geom_smooth(method = "lm")

# Modeling -----

# Using mean surrounding vegetation length as the predictor
tml.msvl.cal.lm <- lm(tot_mat_length_cm ~ mean_surveg_length_cm, data = modsel_cal_plantlvl)

summary(tml.msvl.cal.lm)

### ------------------
### Modeling: L. glaberrima -----

## -------------------
## Do larger plants produce more inflorescences? -----

# Scatterplots -----

ggplot(modsel_gla_plantlvl, aes(y = tot_mat_length_cm, x = n_inflor)) +
        geom_point() + 
        geom_smooth(method = "lm")

# Modeling -----

inflor.totlen.gla.lm <- lm(n_inflor ~ tot_mat_length_cm, data = modsel_gla_plantlvl)

summary(inflor.totlen.gla.lm)

## -------------------
## Do larger plants produce heavier seeds? -----

# Scatterplots -----

ggplot(joined_data_filtered_gla, aes(x = mean_seed_mass, y = total_branch_length)) + 
        geom_point() + 
        geom_smooth(method = "lm")

# Modeling -----

seedmass.totlen.gla.lm <- lm(mean_seed_mass ~ total_branch_length, data = joined_data_filtered_gla)

summary(seedmass.totlen.gla.lm)


## -------------------
## Do taller plants grow in taller vegetation? -----

# Scatterplot -----

# raw heights of both focal plants and surrounding vegetation length 
ggplot(modsel_gla_plantlvl, aes(x = mean_inflor_height_cm, y = mean_surveg_length_cm)) + 
        geom_point() +
        geom_smooth(method = "lm")

# Relative height of focal plant relative to inflorenscence height 
ggplot(modsel_gla_plantlvl, aes(y = rel_inflor_height_to_surveg_height, x = mean_inflor_height_cm)) + 
        geom_point() + 
        geom_smooth(method= "lm")


# Modeling -----

# taller plants tend to occur in taller vegetation (absolute measurements)
mih.msvl.gla.lm <- lm(mean_inflor_height_cm ~ mean_surveg_length_cm, data = modsel_gla_plantlvl)

summary(mih.msvl.gla.lm)

# shorter plants (absolute height) tend to grow in relatively taller surrounding vegetation
relihsvh.msvl.gla.lm <- lm(rel_inflor_height_to_surveg_height ~ mean_inflor_height_cm, data = modsel_gla_plantlvl)

summary(relihsvh.msvl.gla.lm)



## -------------------
## Do larger plants grow in taller vegetation? -----

# Scatterplots -----

ggplot(modsel_gla_plantlvl, aes(x = tot_mat_length_cm, y = mean_surveg_length_cm)) +
        geom_point() + 
        geom_smooth(method = "lm")

ggplot(modsel_gla_plantlvl, aes(x = tot_mat_length_cm, y = surveg_mass_g)) +
        geom_point() + 
        geom_smooth(method = "lm")


# Modeling -----

# Using mean surrounding vegetation length as the predictor
tml.msvl.gla.lm <- lm(tot_mat_length_cm ~ mean_surveg_length_cm, data = modsel_gla_plantlvl)

summary(tml.msvl.gla.lm)



