#### 2019 Field Dispersal Data Analysis 2: Cleaning surrounding vegetation data

# Courtney Van Den Elzen

# Initiated: December 2019
# Latest Update: January 2023 (syntax update for Jan 2023 package version compatibility)


### ------------------
### Load Packages
library(tidyverse)


### ------------------
### Read In Surrounding Vegetation Data 

# surrounding vegetation data
svdf <- read_csv("./Data/Field Study Raw Data - Surrounding Vegetation.csv")


### ------------------
### Clean Surrounding Vegetation Data 

# filter to those that were actually collected (includes those collected that are missing.)
svdf_noomit <- filter(svdf, `Omit from analysis` == "no")

# Change species variable into a factor with three levels: cal, fre, gla
svdf_noomit$Species <- svdf_noomit$Species %>% as.factor() 
levels(svdf_noomit$Species) <- c("cal", "fre", "gla")

# Create ID column
svdf_noomit$Plant_ID <- paste(svdf_noomit$Hub, svdf_noomit$Species, "T", svdf_noomit$Transect, "M", svdf_noomit$Mom, sep = "")

# Make one dataframe that contains the bunch lengths as well as the masses
svdf_mass_len <- dplyr::select(svdf_noomit, "Plant_ID", "Species", "Hub", "Transect", "Mom", "Max length (whole)", "Mass (whole - g)", "Notes")

# clean column names
colnames(svdf_mass_len) <- c("Plant_ID", "Species", "Hub", "Transect", "Mom", "max_length_mm", "Mass_g", "Notes")

# Make another dataframe (transformed to long form) with the individual plant lengths
svdf_ind_len <- dplyr::select(svdf_noomit, "Plant_ID", "Species", "Hub", "Transect", "Mom", "Length (plant 1)", "Length (plant 2)", "Length (plant 3)", "Length (plant 4)", "Length (plant 5)")

# Pivot the data to long format
svdf_ind_len_long <- pivot_longer(svdf_ind_len, c("Length (plant 1)", "Length (plant 2)", "Length (plant 3)", "Length (plant 4)", "Length (plant 5)"), names_to = "Individual Plant Lengths")

# clean column names
colnames(svdf_ind_len_long) <- c("Plant_ID", "Species", "Hub", "Transect", "Mom", "Plant_Num", "Length (mm)")

# Create a df of mean lengths from the 5 separated plants
svdf_len_means <- dplyr::summarize(group_by(svdf_ind_len_long, Plant_ID, Species, Hub, Transect, Mom), mean_len = mean(`Length (mm)`))

# Add the mean lengths back into the original dataframe
svdf_mass_len$mean_length_mm <- svdf_len_means$mean_len

# Calculate surrounding veg density from the mean length, mass, and size of the collection area.
collection_diam <- pi*0.128^2 #area of the circle in mm
cylinder_vol <- collection_diam*(svdf_mass_len$mean_length_mm/10)
svdf_mass_len$veg_density_kg_m2 <- (svdf_mass_len$Mass_g/1000)/cylinder_vol

# Get rid of any NA columns
svdf_mass_len_filt <- svdf_mass_len %>% filter_at(vars(mean_length_mm, Mass_g,max_length_mm,veg_density_kg_m2), any_vars(!is.na(.)))


### ------------------
### Plot Surrounding Vegetation Data 

# Plot the data
ggplot(svdf_mass_len_filt, aes(x=Species, y = mean_length_mm)) + geom_boxplot()
ggplot(svdf_mass_len_filt, aes(x=Species, y = max_length_mm)) + geom_boxplot()
ggplot(svdf_mass_len_filt, aes(x=Species, y = Mass_g)) + geom_boxplot()
ggplot(svdf_mass_len_filt, aes(x=Species, y = veg_density_kg_m2)) + geom_boxplot()


### ------------------
### Write Cleaned Data to File

# Write to a csv for use in later scripts
write_csv(svdf_mass_len_filt, "./Data/Field Study Cleaned Data - Surrounding Veg Mass Max Mean Length.csv")

