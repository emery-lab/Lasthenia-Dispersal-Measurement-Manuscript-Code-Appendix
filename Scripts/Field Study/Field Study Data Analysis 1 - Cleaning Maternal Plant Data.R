#### Field Dispersal Data Analysis 1: Cleaning maternal plant data

# Courtney Van Den Elzen

# Initiated: December 2019
# Latest Update: January 2023 (syntax update for Jan 2023 package version compatibility)


### ------------------
### Load Packages -----
library(tidyverse)


### ------------------
### Read In Maternal Inflorescence Data -----

midf <- read_csv("./Data/Field Study/Field Study Raw Data - Maternal Inflorescences.csv")


### ------------------
### Clean Maternal Inflorescence Data -----

# Get rid of buds (keep only mature inflorescence rows)
midf_nobud <- filter(midf, Size != c("bud")) #also gets rid of NA rows.

# Change species variable into a factor with three levels: cal, fre, gla
midf_nobud$Species <- midf_nobud$Species %>% as.factor() 
levels(midf_nobud$Species) <- c("cal", "fre", "gla")

# Create an ID column
midf_nobud$Plant_ID <- paste(midf_nobud$Hub, midf_nobud$Species, "T", midf_nobud$Transect, "M", midf_nobud$Mom, sep = "")

# Calculate the max number of inflorescences
midf_n_inf <- dplyr::summarize(group_by(midf_nobud, Plant_ID), n_inflor = max(Inflorescence))


### ------------------
### Read In Maternal Plant Data -----

mpdf <- read_csv("./Data/Field Study/Field Study Raw Data - Maternal Plant Traits.csv")

### ------------------
### Clean Maternal Plant Data -----

mpdf <- dplyr::filter(mpdf, `Pins?` == "Yes", !is.na(`Main Stem Length (mm)`))
mpdf$`Pins?` <- NULL

# Subset for pivoting to avoid duplication unnecessary things.
mpdf_sub <- dplyr::select(mpdf, "Full ID", "Main Stem Length (mm)", "Branch 1 Length", "Branch 2 Length", "Branch 3 Length", "Branch 4 Length", "Branch 5 Length", "Branch 5 Length", "Branch 6 Length")

# Pivot into long format
mpdf_long <- pivot_longer(mpdf_sub, 
                          c("Main Stem Length (mm)", "Branch 1 Length", "Branch 2 Length", "Branch 3 Length", "Branch 4 Length", "Branch 5 Length", "Branch 5 Length", "Branch 6 Length"),
                          names_to = "Branch Lengths")

# Change column names
colnames(mpdf_long) <- c("Full ID", "Measurement type", "Length")

# Calculate the total length of each plant
tot_len_df <- dplyr::summarize(group_by(mpdf_long, `Full ID`), Tot_len = sum(Length, na.rm = TRUE))


### ------------------
### Read In Maternal Plant Data (Recast - average flower height is calculable) -----
 
mpdf_r <- read_csv("./Data/Field Study/Field Study Raw Data - Maternal Plant Traits Recast.csv")

# filter out those with no recovered seeds and/or no maternal lengths
mpdf_r <- dplyr::filter(mpdf_r, `Pins?` == "Yes", !is.na(`Full Length Main Branch`))

### ------------------
### Clean Maternal Plant Data (Recast) -----

mpdf_r <- dplyr::filter(mpdf_r, `Pins?` == "Yes", !is.na(`Main Stem Length (mm)`))
mpdf$`Pins?` <- NULL

# filter out some specific columns
mpdf_r_sub <- dplyr::select(mpdf_r, 
                            "Full ID", 
                            "Species", 
                            "Full Length Main Branch", 
                            "Full Length Branch 2", 
                            "Full Length Branch 3", 
                            "Full Length Branch 4", 
                            "Full Length Branch 5", 
                            "Full Length Branch 6",
                            "Full Length Branch 7") %>%
  dplyr::filter(!is.na(`Full Length Main Branch`))


# switch columns and rows (make long form)
mpdf_r_long <- pivot_longer(mpdf_r_sub, 
                          c("Full Length Main Branch", "Full Length Branch 2", 
                            "Full Length Branch 3", "Full Length Branch 4", 
                            "Full Length Branch 5", "Full Length Branch 6",
                            "Full Length Branch 7"),
                          names_to = "Branch Lengths")

# change column names
colnames(mpdf_r_long) <- c("Full ID", "Species", "Measurement type", "Length")

# summary of floral branch length/floral height traits
summary_len_df <- dplyr::summarize(group_by(mpdf_r_long, `Full ID`, Species), 
                                   Var_len = var(Length, na.rm = TRUE), 
                                   Mean_len = mean(Length, na.rm = TRUE),
                                   non_na_count = sum(!is.na(Length)))

# Add total lengths, variance in lengths, and mean length to original dataframe and omit extra columns
mpdf_summary <- dplyr::full_join(mpdf, summary_len_df, by = c("Full ID", "Species")) %>%
     dplyr::select("Full ID", "Species", "Hub", "Transect", "Mom", "Mass (mg)", 
                   "Main Stem Length (mm)", "Mean_len", "Var_len", "non_na_count")


# Add total lengths to original dataframe and omit extra columns
mpdf_totals_means_vars <- 
     full_join(mpdf, tot_len_df, by = "Full ID") %>% 
     full_join(mpdf_summary) %>% dplyr::select("Full ID", "Species", "Hub", 
                                               "Transect", "Mom", "Mass (mg)", 
                                               "Main Stem Length (mm)", "Tot_len", 
                                               "Mean_len", "Var_len", "non_na_count")

### ------------------
### Plot Maternal Plant Data -----

# Plot the data and calculate the correlation between the total length, mean length, and main stem length
ggplot(data = mpdf_totals_means_vars, aes(x = Species, y = Var_len)) +
  geom_boxplot()


### ------------------
### Write Cleaned Data to File -----

# Write to a csv for use in later scripts
#write_csv(mpdf_totals_means_vars, "./Data/Field Study/Field Study Cleaned Data - Maternal Plant Mass Total Mean and Var Lengths.csv")
#write_csv(mpdf_long, "./Data/Field Study/Field Study Cleaned Data - Maternal Plant Branch Lengths Long Form.csv")
#write_csv(mpdf_r_long, "./Data/Field Study/Field Study Cleaned Data - Maternal Plant Full Branch Lengths Recast Long Form.csv")
#write_csv(midf_n_inf, "./Data/Field Study/Field Study Cleaned Data - Maternal Inflorescences.csv")

