#### Drop Tube Data Analysis 2: Data Analysis and Visualization

# Courtney Van Den Elzen

# Original: 2019

# Latest Update: January 2023 (syntax update for Jan 2023 version compatibility)

# Description: This script reformats the dataframes to being mean velocity 
#              instead of having all velocity measurements. It then builds plots 
#              and linear models of seed mass and species effects on velocity 
#              (both raw measurements and means)


### ------------------
### Load Packages -------

library(lme4)
library(tidyverse)
library(nlme)


### ------------------
### Read In Data  -------

dt_data <- read.csv("./Data/Drop Tube Study/Drop Tube Cleaned Data - Expanded.csv")

### ------------------
### Clean Data -------

# drop NAs from the dataframe
dt_data <- drop_na(dt_data) 

# Create data frame of the average velocity value for each video

# In "Drop Tube Data Analysis 1.R" the first and last three measurements were 
# dropped due to boundary issues. This is the average of the measurements that 
# are left after that boundary cropping
dt_data_means <- dplyr::summarize(
  group_by(dt_data, 
           Video_Code, 
           Species, 
           Plant_ID, 
           Seed_Type, 
           Seed_Mass_mg), 
  mean_velocity = mean(Velocity_Vals_m_sec))

# Break down means df by species 
dt_data_means_cal <- filter(dt_data_means, Species == "cal")
dt_data_means_fre <- filter(dt_data_means, Species == "fre")
dt_data_means_gla <- filter(dt_data_means, Species == "gla")

# Further break down by seed type as well
dt_data_means_cal_disc <- filter(dt_data_means, Species == "cal", Seed_Type == "d")
dt_data_means_fre_disc <- filter(dt_data_means, Species == "fre", Seed_Type == "d")

### ------------------
### Data Visualization -------

# Species by mean velocity boxplot
ggplot(dt_data_means, aes(Species, mean_velocity)) + geom_boxplot()

# Seed mass by species boxplot
ggplot(dt_data_means, aes(Species, Seed_Mass_mg)) + geom_boxplot()

# Seed mass by seed type boxplot by species
ggplot(dt_data_means, aes(Seed_Type, Seed_Mass_mg)) + geom_boxplot() + facet_wrap(~Species)

# Mean velocity by species boxplot
ggplot(dt_data_means, aes(mean_velocity)) + geom_histogram() + facet_wrap(~Species)

# Seed mass by species boxplot
ggplot(dt_data_means, aes(Seed_Mass_mg)) + geom_histogram() + facet_wrap(~Species)

# Scatterplot of velcity by seed mass for each species 
ggplot(dt_data_means, aes((log(Seed_Mass_mg)), log(mean_velocity), color=Species)) + 
  geom_point() +
  scale_color_manual(labels = c("L. californica", "L. fremontii", "L. glaberrima"), values=c("#999999", "#E69F00", "#56B4E9")) +
  theme_classic() +
  geom_smooth(method='lm') +
  xlab("log(Seed Mass (mg))") + 
  ylab("log(Mean Drop Velocity (m/s))") +
  theme(axis.title.x = element_text(size=14),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        axis.title.y = element_text(size=14),
        legend.title = element_text(face = "italic", size=12), 
        legend.text = element_text(face = "italic", size=12), 
        legend.position = c(0.82, 0.25), 
        legend.background = element_rect(linetype="solid", colour ="black"))

# Velocity by species boxplot for all species
ggplot(dt_data_means, aes(Species, mean_velocity)) + 
  geom_boxplot() + 
  labs(y = "Mean Drop Velocity (m/s)") +
  #scale_y_continuous(breaks=seq(0,59,10)) + 
  scale_x_discrete(labels=c("L. californica", "L. fremontii", "L. glaberrima")) +
  theme_classic() +
  theme(axis.text.x = element_text(face = "italic"))

# Velocity by seed mass scatterplot for L. californica disc seeds
ggplot(dt_data_means_cal_disc, aes(Seed_Mass_mg, mean_velocity)) + 
  geom_point() +
  geom_smooth(method = "lm", se = T) +
  xlab("Seed Mass (mg)") + 
  ylab("Mean Drop Velocity (m/s)") +
  ggtitle("Lasthenia californica") +
  theme(plot.title = element_text(color = "red", face="italic"))

# Velocity by seed mass scatterplot for L. fremontii disc seeds
ggplot(dt_data_means_fre_disc, aes(Seed_Mass_mg, mean_velocity)) +
  geom_point() +
  geom_smooth(method = "lm", se = T) +
  xlab("Seed Mass (mg)") + 
  ylab("Mean Drop Velocity (m/s)") +
  ggtitle("Lasthenia fremontii") +
  theme(plot.title = element_text(color = "chartreuse4", face="italic"))

# Velocity by seed mass scatterplot for L. glaberrima "disc" seeds
ggplot(dt_data_means_gla, aes(Seed_Mass_mg, mean_velocity)) +
  geom_point() +
  geom_smooth(method = "lm", se = T) +
  xlab("Seed Mass (mg)") + 
  ylab("Mean Drop Velocity (m/s)") +
  ggtitle("Lasthenia glaberrima") +
  theme(plot.title = element_text(color = "blue", face="italic"))

### ------------------
### Statistical Modeling -------

# Fit linear models for each of the species independently on the mean velocity values

# L. californica simple linear model
modmeans.cal <- lm(mean_velocity ~ Seed_Mass_mg, data = dt_data_means_cal_disc)

summary(modmeans.cal)
par(mfrow = c(2, 2))
plot(modmeans.cal)

# L. fremontii simple linear model
modmeans.fre <- lm(mean_velocity ~ Seed_Mass_mg, data = dt_data_means_fre_disc)

summary(modmeans.fre)
par(mfrow = c(2, 2))
plot(modmeans.fre)

# L. glaberrima simple linear model
modmeans.gla <- lm(mean_velocity ~ Seed_Mass_mg, data = dt_data_means_gla)
 
summary(modmeans.gla)
par(mfrow = c(2, 2))
plot(modmeans.gla)

# Drop the outlier and replot the data
dt_data_means_gla_noOL <- subset(dt_data_means_gla, Seed_Mass_mg < 0.4)
 
ggplot(dt_data_means_gla_noOL, aes(Seed_Mass_mg, mean_velocity)) +
  geom_point() +
  geom_smooth(method = "lm", se = T) +
  xlab("Seed Mass (mg)") + 
  ylab("Mean Drop Velocity (m/s)") +
  ggtitle("Lasthenia glaberrima") +
  theme(plot.title = element_text(face="italic"))

