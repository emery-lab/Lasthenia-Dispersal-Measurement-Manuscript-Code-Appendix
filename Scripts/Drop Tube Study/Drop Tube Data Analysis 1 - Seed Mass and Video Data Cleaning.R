#### Drop Tube Data Analysis 1: Data Wrangling and Categorical Set Up

# Courtney Van Den Elzen

# Original: 2019

# Latest Update: January 2023 (syntax update for Jan 2023 version compatibility)

# Description: Read in and clean the drop tube seed mass and video velocity data
#              and metadata

### ------------------
### Load Packages -------

library(tidyverse)
library(knitr)
library(rsalad)


### ------------------
### Read In Data  -------

# Import metadata and seed data
metad <- read.csv("./Data/Drop Tube Study/Drop Tube Raw Data - Meta Data.csv")
massd <- read.csv("./Data/Drop Tube Study/Drop Tube Raw Data - Seed Data.csv")

### ------------------
### Clean Data -------

# Wrangle the data into one data frame

# Check if the two datasheets line up - they do
all.equal(metad$Video_Code, massd$Video_Code)

# Merge the data frames by the video code
meta_mass_merge <- merge(metad, 
                         massd, 
                         by="Video_Code", 
                         all=F)

# Create a new column that is a plant-level ID
meta_mass_merge$Plant_ID <- paste(meta_mass_merge$Set, 
                                  meta_mass_merge$Species,
                                  meta_mass_merge$Number, 
                                  "m", 
                                  meta_mass_merge$Mom, 
                                  "n", 
                                  meta_mass_merge$Neighbour, 
                                  sep='')


# Extract just important columns - all ID info and mass measurements.
metad_base <- meta_mass_merge[,c(1:8,18,20)]

# Extract the pixels per mm
pix_mm <- metad$Pixels_mm

# Extract number of different seeds in the data frame
nseeds <- dim(metad_base)[1]


### ------------------
### Read In Velocity Data Files -------

# initialize a vector to hold data vector length values
datlen2 <- c()
all_data_vec <- c()
data_index_vec <- c()

# Get max length of a velocity vector
for (i in 1:nrow(metad_base)){
  
  # import data from given video
  vidname <- paste("./Data/Drop Tube Study/Velocity Data Files/", 
                   metad_base$Video_Code[i], 
                   "_Data.txt", 
                   sep='')
  
  if (i == 10){
    print(vidname)
  }
  
  # read in data from video (vidname) - units weird from MATLAB conversion mistake
  data <- as.vector(t(read.csv(vidname, header=F)))
  datlen <- length(data)
  
  if (i == 10){
    print("data - weird units (fixed further down)")
    print(data)
    print("datlen - length of the velocity vector")
    print(datlen)
  }
  
  #' take out first and last three data points because of boundary issues
  data2 <- data[-c(1, 2, 3, (datlen-2), (datlen-1), datlen)]
  data2len <- length(data2)
  
  # Correct the unit mistake made in MATLAB (wrong conversion calculation - easier to fix here)
  # Mistake made in MATLAB was multiplying by pix/mm instead of dividing, giving units of 
  # pix^2 / sec*mm instead of mm/sec. This corrects units to mm/sec
  data2_corr <- data2 / pix_mm[i] / pix_mm[i]
  
  if (i == 10){
    print("data2_corr - velocity measurements units fixed")
    print(data2_corr)
  }
  
  # Plop the data into the data3 vector
  data3 <- c()
  data3[1:61] <- NA # 61 is the longest vector length for velocity values - should not be hard coded
  data3[1:data2len] <- data2_corr
  
  # calculate new vector length, add to vector of length values
  datlen2 <- c(datlen2, length(data2))
  
  if (i == 10){
    print("data2_corr - velocity measurements units fixed")
    print(data2_corr)
  }
  
  # Make a long vector of all of the cleaned velocity data points. The units on this are mm/sec
  all_data_vec <- c(all_data_vec, data3)
  
  
  # Make an indexing vector with 1-61 for each video (to house the max = 61 velocity measurements kept.)
  data_index_vec <- c(data_index_vec, seq(from=1, to=61))
  
}


### ------------------
### Clean Velocity Data -------

# New column with the lengths of the data vectors from each video
metad_base$Num_Vals <- datlen2

# Expand the data frame for one row for each entry from each data vector
metad_base.expanded <- metad_base[rep(row.names(metad_base), each=61), 1:10]

# Add new column with the velocity data
rownames(metad_base.expanded) <- NULL
metad_base.expanded$Velocity_Vals_m_sec <- all_data_vec/1000
metad_base.expanded$Data_Index <- data_index_vec

# Get rid of NAs
metad_base.expanded.naomit <- na.omit(metad_base.expanded)

# Summarize the expanded data frame by mean velocity value

# Create the means vector
metad_base.exp.means <- na.omit(metad_base.expanded) %>%
  dplyr::group_by(Video_Code) %>%
  dplyr::summarise(Velocity_Vals_m_sec = mean(Velocity_Vals_m_sec))

# Merge with the original data frame to add in all other variables
metad_base.exp.means <- merge(metad_base, 
                              metad_base.exp.means, 
                              by = "Video_Code")


#write.csv(metad_base.expanded, "./Data/Drop Tube Study/Drop Tube Cleaned Data - Expanded.csv")
#write.csv(metad_base.exp.means, "./Data/Drop Tube Study/Drop Tube Cleaned Data - Expanded Means.csv")
