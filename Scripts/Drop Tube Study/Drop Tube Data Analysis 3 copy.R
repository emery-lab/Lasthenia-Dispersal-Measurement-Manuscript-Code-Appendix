#' # Drop Tube Data Analysis Part 3
#' 
#' ### Cleaning continuous trait data
#' 
#' #### Courtney Van Den Elzen 
#'
#' #### 2018-2019 Study
#' 

#' ### Libraries
#+ load-libs, include = TRUE
library(lme4)
library(tidyverse)
library(nlme)
library(cubature)

#' Load in the original data (all velocity measurements)
dt_data <- read_csv("/Users/Courtney/Documents/Thesis Chapter 2/Drop Tube Study/Data/DropTubeAllDataExpanded.csv")
dt_data$X1 <- NULL
dt_data <- dt_data %>% mutate_at(vars(Number, Mom, Neighbour, Seed_Number), factor) 

#' Load in the means data (average velocity measurements)
dt_data_means <- read_csv("/Users/Courtney/Documents/Thesis Chapter 2/Drop Tube Study/Data/DropTubeAllDataExpanded_Means.csv")
dt_data_means$X1 <- NULL
dt_data_means <- dt_data_means %>% mutate_at(vars(Number, Mom, Neighbour, Seed_Number), factor) 
dt_data_means %>% dplyr::group_by(Plant_ID, Seed_Type) %>% dplyr::summarize(n()) %>% View
#' Load in the seed trait data
traits <- read_csv("/Users/Courtney/Documents/Thesis Chapter 2/Drop Tube Study/Data/Drop Tube Seed Trait Data.csv")

#' Correct species labels that are incorrect
traits <- traits %>% mutate(Species = case_when(Other_species == "no" ~ Species,
                                                Other_species == "cal" ~ "L. californica",
                                                Other_species == "fre" ~ "L. fremontii",
                                                Other_species == "gla" ~ "L. glaberrima"))


#' Filter out seeds with trait issues
traits <- traits %>% filter(Angle_problem == "no" & Other_issue == "no")

#' Get a list of all of the video codes in the trait dataset
trait_VC <- as.data.frame(unique(dt_data$Video_Code))
colnames(trait_VC) <- "Video_Code"

#' This is just subsetting the dt_data so the video codes exactly match those from the trait data
#' Some seed images had issues that led them to being omitted from the trait dataset. The original trait data
#' was only collected on seeds that I had successfully gotten velocity data for so the dt data is a subset 
#' of the pruned trait data
dt_data_means <- dt_data_means %>% merge(trait_VC, by = "Video_Code") %>% droplevels()

#' Subset to the traits that are important for the downstream analysis
traits_sub <- traits[,c("Video_Code", "Species", "Orientation", "Pappus_base_width", "Pappus_tip_width", "Achene_tip_width", "Achene_mid_width", "Achene_length_lin", "Achene_twopt_length", "Achene_twopt_angle", "Achene_max_width", "Pappus_perp_max_length")]

#' Split into front and side traits
traits_sub_f <- arrange(filter(traits_sub, Orientation == "Front"), Video_Code) %>% mutate(Orientation = NULL)
traits_sub_s <- arrange(filter(traits_sub, Orientation == "Side"), Video_Code) %>% mutate(Orientation = NULL)

#' Rename the columns so they're unique to front and side
colnames(traits_sub_f) <- c("Video_Code", "Species", "Pappus_base_width_f", "Pappus_tip_width_f", "Achene_tip_width_f", "Achene_mid_width_f", "Achene_length_lin_f", "Achene_twopt_length_f", "Achene_twopt_angle_f", "Achene_max_width_f", "Pappus_perp_max_length_f")
colnames(traits_sub_s) <- c("Video_Code", "Species", "Pappus_base_width_s", "Pappus_tip_width_s", "Achene_tip_width_s", "Achene_mid_width_s", "Achene_length_lin_s", "Achene_twopt_length_s", "Achene_twopt_angle_s", "Achene_max_width_s", "Pappus_perp_max_length_s")

traits_merge <- merge(traits_sub_f, traits_sub_s, by=c("Video_Code","Species"))
traits_merge[1:3,]

#' Traits to add: "Achene.twopt.angle.dev.max", "pap.tip.cyp.max.area.ratio", "pap.tip.pap.base.area.ratio", "Pap.len.ach.len.ratio", "pap.tip.ecc", "coneSA"
#' 
#' Add in some modified base traits for construction of other traits
#' 
#' **Average Achene Length**
#+
Ach.len.lin.mean <- apply(cbind(traits_sub_f$Achene_length_lin_f, traits_sub_s$Achene_length_lin_s), 1, mean)

#' **Average Pappus Length**
#+
Pap.perp.max.len.mean <- apply(cbind(traits_sub_f$Pappus_perp_max_length_f, traits_sub_s$Pappus_perp_max_length_s), 1, mean)

#' **Two-Point Angle Deviation From 180**
#' 
Achene_twopt_angle_dev_f <- 180-traits_sub_f$Achene_twopt_angle_f
Achene_twopt_angle_dev_s <- 180-traits_sub_s$Achene_twopt_angle_s

#' **Max Angle**
#+
Achene_twopt_angle_dev_max <- pmax(Achene_twopt_angle_dev_f, Achene_twopt_angle_dev_s)
Achene_twopt_angle_dev_min <- pmin(Achene_twopt_angle_dev_f, Achene_twopt_angle_dev_s)

#' **Pappus Length Ratios**
#+
Pap.len.ach.len.ratio <- Pap.perp.max.len.mean/Ach.len.lin.mean

#' Defining pi for area calculations
#+
pi = 3.14159265359

#' **Area Traits**
#+
pappus_tip_area <- pi*(traits_sub_f$Pappus_tip_width_f/2)*(traits_sub_s$Pappus_tip_width_s/2)
pappus_base_area <- pi*(traits_sub_f$Pappus_base_width_f/2)*(traits_sub_s$Pappus_base_width_s/2)
cypsela_max_area <- pi*(traits_sub_f$Achene_max_width_f/2)*(traits_sub_s$Achene_max_width_s/2)

#' **Area Ratio Traits**
#+
pap.tip.cyp.max.area.ratio <- pappus_tip_area/cypsela_max_area
pap.tip.pap.base.area.ratio <- pappus_tip_area/pappus_base_area

#' **Pappus Tip Eccentricity**
#' 
#' Find the semi-major and minor axes for each pappus tip ellipse
#+
major_ax_tip <- apply(cbind(traits_sub_f$Pappus_tip_width_f, traits_sub_s$Pappus_tip_width_s), 1, max)/2
minor_ax_tip <- apply(cbind(traits_sub_f$Pappus_tip_width_f, traits_sub_s$Pappus_tip_width_s), 1, min)/2

#' **Pappus Base Eccentricity**
#' 
#' Find the semi-major and minor axes for each pappus base ellipse (this is the "radius" in a given direction) 
#+
major_ax_base <- apply(cbind(traits_sub_f$Pappus_base_width_f, traits_sub_s$Pappus_base_width_s), 1, max)/2
minor_ax_base <- apply(cbind(traits_sub_f$Pappus_base_width_f, traits_sub_s$Pappus_base_width_s), 1, min)/2

#' Calculate eccentricity
#+
pap_tip_ecc <- sqrt(1-((minor_ax_tip^2)/(major_ax_tip^2)))

#' **Cone Surface Area**
#' 
#' use the defined elliptic cone formulae here: http://mathworld.wolfram.com/EllipticCone.html and then modified the code above.
#+
coneA <- function(a,b,h) {
  require(cubature)
  integrand <- function(x,a,b){
    u <- x[1]
    v <- x[2]
    E <- (h^2 + (a*cos(v))^2 + (b*sin(v))^2)/(h^2)
    F1 <- (a^2 - b^2)*(h-u)*(cos(v)*sin(v))/(h^2)
    G <- ((h-u)^2)*((a*sin(v))^2 + (b*cos(v))^2)/(h^2)
    return(sqrt(E*G-F1^2))
  }
  adaptIntegrate(integrand, lowerLimit=c(0,0), upperLimit=c(h,2*pi), a=a, b=b)$integral
}

ach.cone.data <- data.frame("a" = major_ax_base, "b" = minor_ax_base, "h" = Ach.len.lin.mean)

#' Make the data vector
#+
coneSA <- mapply(coneA, a = ach.cone.data$a, b = ach.cone.data$b, h = ach.cone.data$h)

#' Add all resultant traits to the dataframe
traits_merge$pappus_tip_area <- pappus_tip_area
traits_merge$Achene.twopt.angle.dev.max <- Achene_twopt_angle_dev_max 
traits_merge$pap.tip.cyp.max.area.ratio <-  pap.tip.cyp.max.area.ratio 
traits_merge$pap.tip.pap.base.area.ratio <- pap.tip.pap.base.area.ratio
traits_merge$Pap.len.ach.len.ratio <- Pap.len.ach.len.ratio
traits_merge$pap.tip.ecc <- pap_tip_ecc 
traits_merge$coneSA <- coneSA

#' Merge the seed metadata and velocity data with the trait data
dt_traits_merge <- merge(dt_data_means, traits_merge, by="Video_Code") %>% mutate(Species.x = NULL) %>% rename(Species = Species.y)


#knitr::spin("/Users/Courtney/Documents/Thesis Chapter 2/Drop Tube Study/Analysis/Drop Tube Data Analysis 3.R")
