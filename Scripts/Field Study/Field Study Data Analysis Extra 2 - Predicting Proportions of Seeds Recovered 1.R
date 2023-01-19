#### 2019 Field Dispersal Data Analysis Extra 2: Predicting Proportions of Seed Recovered 1

# Courtney Van Den Elzen

# Initiated: October 2019
# Latest Update: January 2023 (syntax update for Jan 2023 package version compatibility)

# Description: Clean and combine the 2018 and 2019 field neighbour datasets 
#              (with seed counts by diameter) to train a predictive model


### ------------------
### Load Packages -----

library(tidyverse) # for cleaning, plotting
library(MASS) # extra stats stuff
library(doParallel) # to do things in parallel... not sure if this is actually used.


### ------------------
### Read In 2018 (pilot study) Neighbor Data -----

# Load data
nseed_data_18 <- read_csv("./Data/Field Study/Field Study Raw Data - 2018 Neighbour Data.csv") 

### ------------------
### Clean 2018 Neighbor Data -----

colnames(nseed_data_18)[1] <- "Hub"
nseed_data_18$Mat_ID <- as.factor(paste(nseed_data_18$Hub, nseed_data_18$Number, nseed_data_18$Species, "t1", "m", nseed_data_18$Mom, sep = ""))
nseed_data_18$Indiv_ID <- as.factor(paste(nseed_data_18$Mat_ID,"n", nseed_data_18$Neighbour, sep = ""))


# Create ID and year columns
nseed_data_18$Year <- as.factor("2018")
nseed_data_18$MID_Year <- as.factor(paste(nseed_data_18$Year, nseed_data_18$Mat_ID, sep = ""))
nseed_data_18$IID_Year <- as.factor(paste(nseed_data_18$Year, nseed_data_18$Indiv_ID, sep = ""))

# Filter out NA rows for diameter measurements and problematic rows
nseed_data_18_1flr <- droplevels(dplyr::filter(nseed_data_18, 
                                               !is.na(Diameter) & 
                                                 Omit_from_analysis == "no"))

# Subset to the necessary columns
nseed_data_18_1flr <- 
  nseed_data_18_1flr %>% 
  dplyr::select(Species, Mat_ID, Indiv_ID, Year, MID_Year, IID_Year, Diameter, 
                Viable_Disc_Seed_Count, Viable_Ray_Seed_Count)

#' Fix column names
colnames(nseed_data_18_1flr) <- c("Species", "Mat_ID", "Indiv_ID", "Year", 
                                  "MID_Year", "IID_Year", "Diameter", 
                                  "Disc_tube_count", "Ray_tube_count")


### ------------------
### Read In 2019 (pilot study) Neighbor Data -----

# Load the data 
nseed_data_19 <- read_csv("./Data/Field Study/Field Study Raw Data - Neighbour Seed Counts.csv") %>%
  dplyr::mutate(species = case_when(species == "L. californica" ~ "cal",
                                    species == "L. fremontii" ~ "fre",
                                    species == "L. glaberrima" ~ "gla"))

colnames(nseed_data_19)


### ------------------
### Clean 2019 Neighbor Data -----

# rename columns & levels of species to match 2018 data

colnames(nseed_data_19) <- c("Species", "Hub", "Transect", "Mom", "Neighbour", "Plant_ID", "Count_type", 
  "Viable_Ray_Seed_Count", "Viable_Disc_Seed_Count", "Total_seeds", "Ray_tube_count", 
  "Disc_tube_count", "Flowers_in_count", "Diameter", "Perfect_separation", "Notes")

nseed_data_19

# Make Year and ID columns
nseed_data_19$Mat_ID <- as.factor(paste(nseed_data_19$Hub, 1, nseed_data_19$Species, "t", nseed_data_19$Transect, "m", nseed_data_19$Mom, sep = ""))
nseed_data_19$Indiv_ID <- as.factor(paste(nseed_data_19$Hub, 1, nseed_data_19$Species, "t", nseed_data_19$Transect, "m", nseed_data_19$Mom, "n", nseed_data_19$Neighbour, sep = ""))
nseed_data_19$Year <- as.factor("2019")
nseed_data_19$MID_Year <- as.factor(paste(nseed_data_19$Year, nseed_data_19$Mat_ID))
nseed_data_19$IID_Year <- as.factor(paste(nseed_data_19$Year, nseed_data_19$Indiv_ID, sep = ""))



# Filter to those individuals with diameter measurements
nseed_data_19_diam <- 
  nseed_data_19 %>% 
  dplyr::filter(!is.na(Diameter)) %>% 
  dplyr::filter(Perfect_separation == "Yes")

nseed_data_19_diam

# Subset to useful columns
nseed_data_19_diam <- 
  nseed_data_19_diam %>% 
  dplyr::select("Species", "Mat_ID", "Indiv_ID", "Year", "MID_Year", "IID_Year", 
                "Diameter", "Disc_tube_count", "Ray_tube_count")


# Combine the 2018 and 2019 data and plot it
nseed_data_all_good <- bind_rows(nseed_data_18_1flr, nseed_data_19_diam)

# Make a column of seed totals
nseed_data_all_good <- 
  nseed_data_all_good %>% 
  dplyr::mutate(Total_count = Disc_tube_count + Ray_tube_count)

### ------------------
### Plot the Data -----

# boxplot of seed counts by MID year
ggplot(nseed_data_all_good, aes(x = MID_Year, y = Disc_tube_count, colour = Species)) +
  geom_boxplot()


# Plot of the simple linear between the count of viable disc seeds and the inflorescence diameter
ggplot(nseed_data_all_good, aes(x=Diameter, y=Disc_tube_count, colour = Species)) +
  geom_jitter() +
  geom_smooth(method=lm, se=T)+
  theme_classic() +
  theme(legend.position="bottom", 
        legend.text = element_text(face = "italic"), 
        legend.background = element_rect(colour ="grey25")) +
  labs(y = "Disc tube count", x = "Inflorescence diameter (cm)") +
  scale_x_continuous(limits=c(0, 11), breaks=seq(0,11,2))

### ------------------
### Write Combined 2018/2019 Data to File -----

write_csv(nseed_data_all_good, file = "./Data/Field Study/Field Study Cleaned Data - 2018 and 2019 Neighbour Data Cleaned and Combined.csv")


