#### 2019 Field Dispersal Data Analysis 3: Cleaning Windspeed Data

# Courtney Van Den Elzen

# Initiated: December 2019
# Latest Update: January 2023 (syntax update for Jan 2023 package version compatibility)


### ------------------
### Load Packages
library(tidyverse)


### ------------------
### Read In Wind Speed Data

hwdf <- read_csv("./Data/Field Study Raw Data - Wind Speed.csv")
hwdf$Notes <- NULL
hwdf <- hwdf %>% filter_all(any_vars(!is.na(.)))
colnames(hwdf) <- c("Full ID", "Species", "SppID", "Hub", "Transect", "Mom", "Wind_speed_1m", "Wind_speed_10cm", "Time", "Day")

# Proportional wind speed
hwdf$prop_speed <- hwdf$Wind_speed_10cm/hwdf$Wind_speed_1m


### ------------------
### Clean Wind Speed Data

# calculate the mean wind speeds at 1m, 10cm, and the mean proportional speed
hwdf_means <- summarize(group_by(hwdf, `Full ID`, SppID, Hub, Transect, Mom, Time, Day), 
                        mean_speed_1m = mean(Wind_speed_1m), 
                        mean_speed_10cm = mean(Wind_speed_10cm),
                        mean_prop_speed = mean(prop_speed))


### ------------------
### Plot Surrounding Vegetation Data 

ggplot(hwdf, aes(x = SppID, y = Wind_speed_1m)) + geom_boxplot()
ggplot(hwdf, aes(x = SppID, y = Wind_speed_10cm)) + geom_boxplot()
ggplot(hwdf, aes(x = SppID, y = prop_speed)) + geom_boxplot()


### ------------------
### Write Cleaned Data to File

# Write to a csv for use in later scripts
write_csv(hwdf, "./Data/Field Study Cleaned Data - Wind Speed.csv")
write_csv(hwdf_means, "./Data/Field Study Cleaned Data - Wind Speed Means.csv")

