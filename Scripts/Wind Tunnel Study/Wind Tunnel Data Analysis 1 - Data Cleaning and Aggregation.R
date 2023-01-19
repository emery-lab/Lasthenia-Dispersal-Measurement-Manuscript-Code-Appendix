#### Wind Tunnel Data Analysis 1: Trait Data Cleaning and Aggregation

# Courtney Van Den Elzen

# Original: 2018 

# Latest Update: January 2023 (syntax update for Jan 2023 version compatibility)

# Description: Cleaning & aggregating data and calculating derived traits

# Guide to shorthands used in naming of traits etc.:
# - pap = pappus
# - dist = distance
# - ach = achene
# - len = length
# - lin = linear
# - fnt = front
# - twopt = two-point
# - ecc = eccentricity
# - ax = axis

### ------------------
### Load Packages -------

library(cubature)
library(lme4)
library(tidyverse)


### ------------------
### Read In Data -------
 
# Read in the distance data for all found seeds
dists <- read.csv("./Data/Wind Tunnel Study/Wind Tunnel Raw Data - Distance and Categorical Data Found Seeds.csv")

# Read in the trait data for all found seeds
traits_mod <- read.csv("./Data/Wind Tunnel Study/Wind Tunnel Raw Data - Seed Trait Data Found Seeds.csv")

### ------------------
### Clean Data Part 1 -------

## -------------
## Omitting Data Points (with explanations) -------

# Subset to the main study categorical and distance data 

# Seed G40S1 is likely gracilis so we're removing it 
dists_main <- subset(dists, Study != "side", Tube.ID != c("G40S1"))


# Subset out all images that have angle issues or other issues
traits_mod2 <- 
  traits_mod %>% 
  dplyr::filter(Problem.with.angle. == "no" & Other.issue. == "no")


# Get the set of unique Tube IDs from the trait data
traits_TubeID <- as.data.frame(unique(traits_mod$Tube.ID))

# Change column name for merging 
colnames(traits_TubeID) <- "Tube.ID"

dists2_main <- 
  dists_main %>% 
  merge(traits_TubeID, by = "Tube.ID") %>% 
  droplevels()

## -------------
## Data Transformations -------

# based on transformations of other trait values 

# Scale and center seed mass data
dists2_main$Seed_Weight_std <- as.vector(scale(dists2_main$Seed_Weight))

# Create "hypoteneuse" distance variable, which is just using pythagoras on the 
# x and y distances
dists2_main$Total.Hypoteneuse.Distance <- 
  sqrt(dists2_main$Total_Y_Distance^2 + 
         dists2_main$Total_X_Distance_zero_mean^2)


### ------------------
### Clean Data Part 2 -------

## -------------
## Data Transformations -------

# This section contains regular data cleaning but also subsetting of data into 
# "Front" and "Side" measurements for future calculation of derived traits. 
# Here, Front and Side are the side of the seed imaged (~90 degrees apart)
 
traits_mod$Who <- NULL
traits_mod$Redo <- NULL
traits_mod$Date <- NULL

# Break trait data into front and side measurements
traits_mod_f <- arrange(filter(traits_mod, Orientation == "Front"), Tube.ID)
traits_mod_s <- arrange(filter(traits_mod, Orientation == "Side"), Tube.ID)

# Create data frame to store aggregate measures, derived measures, and averages 
# between fronts and sides.
 
# Grabs all of the relevant columns from the data frame of FRONT traits
traits_mod_f_agg <- 
  traits_mod_f[,c("Tube.ID", "Species", "Wild.Grnhs", "Orientation", "Pappus_base_width", 
                  "Pappus_tip_width", "Achene_tip_width", "Achene_mid_width", 
                  "Achene_length_lin", "Achene_twopt_length", "Achene_twopt_angle", 
                  "Achene_max_width_.eyeballed.", "Number_awns", "Number_scale_tips", 
                  "Pappus_perp_max_length", "Scales_max_length", "Scales_min_length", 
                  "Awns_max_length", "Awns_min_length")]

# Grabs all of the relevant columns from the data frame of SIDE traits
traits_mod_s_agg <- 
  traits_mod_s[,c("Tube.ID", "Species", "Wild.Grnhs", "Orientation", 
                  "Pappus_base_width", "Pappus_tip_width", "Achene_tip_width", 
                  "Achene_mid_width", "Achene_length_lin", "Achene_twopt_length", 
                  "Achene_twopt_angle", "Achene_max_width_.eyeballed.", "Number_awns", 
                  "Number_scale_tips", "Pappus_perp_max_length", "Scales_max_length", 
                  "Scales_min_length", "Awns_max_length", "Awns_min_length")]

# Merge based on the seed ID
traits_mod_agg <- merge(traits_mod_f_agg, traits_mod_s_agg, by=c("Tube.ID","Species",
                                                                 "Wild.Grnhs"), all=TRUE)


# Clean up the column names

colnames(traits_mod_agg) <- c("Tube.ID", "Species", "Wild.Grnhs", "Orientation_f", 
                              "Pappus_base_width_f", "Pappus_tip_width_f", 
                              "Achene_tip_width_f", "Achene_mid_width_f", 
                              "Achene_length_lin_f", "Achene_twopt_length_f", 
                              "Achene_twopt_angle_f", "Achene_max_width_f", 
                              "Number_awns_f", "Number_scale_tips_f", 
                              "Pappus_perp_max_length_f", "Scales_max_length_f", 
                              "Scales_min_length_f", "Awns_max_length_f" , 
                              "Awns_min_length_f", "Orientation_s","Pappus_base_width_s", 
                              "Pappus_tip_width_s", "Achene_tip_width_s", "Achene_mid_width_s", 
                              "Achene_length_lin_s", "Achene_twopt_length_s", "Achene_twopt_angle_s", 
                              "Achene_max_width_s", "Number_awns_s", "Number_scale_tips_s", 
                              "Pappus_perp_max_length_s", "Scales_max_length_s", 
                              "Scales_min_length_s", "Awns_max_length_s", "Awns_min_length_s")

# Add in mass data (standardized and not)
traits_mod_agg$Seed.weight.std <- dists2_main$Seed_Weight_std
traits_mod_agg$Seed.weight <- dists2_main$Seed_Weight

# Add in pappus hair count data (awn and scale number data)
traits_mod_agg$Awn.count <- pmax(traits_mod_agg$Number_awns_f, traits_mod_agg$Number_awns_s)
traits_mod_agg$Scale.count <- pmax(traits_mod_agg$Number_scale_tips_f, traits_mod_agg$Number_scale_tips_s)

### ------------------
### Clean Data Part 3 -------

## -------------
## Min and Max of Each Width Value -------

# Pappus tip width max and min
traits_mod_agg$Pap.tip.width.max <- pmax(traits_mod_agg$Pappus_tip_width_f, 
                                         traits_mod_agg$Pappus_tip_width_s)
traits_mod_agg$Pap.tip.width.min <- pmin(traits_mod_agg$Pappus_tip_width_f, 
                                         traits_mod_agg$Pappus_tip_width_s)

# Pappus base width max and min
traits_mod_agg$Pap.base.width.max <- pmax(traits_mod_agg$Pappus_base_width_f, 
                                          traits_mod_agg$Pappus_base_width_s)
traits_mod_agg$Pap.base.width.min <- pmin(traits_mod_agg$Pappus_base_width_f, 
                                          traits_mod_agg$Pappus_base_width_s)

# Achene tip width max and min
traits_mod_agg$Ach.tip.width.max <- pmax(traits_mod_agg$Achene_tip_width_f, 
                                         traits_mod_agg$Achene_tip_width_s)
traits_mod_agg$Ach.tip.width.min <- pmin(traits_mod_agg$Achene_tip_width_f, 
                                         traits_mod_agg$Achene_tip_width_s)


# Achene max width max
traits_mod_agg$Ach.max.width.max <- pmax(traits_mod_agg$Achene_max_width_f, 
                                         traits_mod_agg$Achene_max_width_s)


## -------------
## Pappus Tip Width Ratios -------

# Creating ratios of pappus base and tip width (one ratio per "direction", 
# ratio>1 == pappus tip > pappus base)
Pap.tip.base.width.ratio_f <- traits_mod_f$Pappus_tip_width/traits_mod_f$Pappus_base_width
Pap.tip.base.width.ratio_s <- traits_mod_s$Pappus_tip_width/traits_mod_s$Pappus_base_width

# Finding max and min of the two calculated ratios
Pap.tip.base.width.ratio.max <- pmax(Pap.tip.base.width.ratio_f, Pap.tip.base.width.ratio_s)
Pap.tip.base.width.ratio.min <- pmin(Pap.tip.base.width.ratio_f, Pap.tip.base.width.ratio_s)

# Ratio of pappus tip to achene max width (one ratio per "direction", 
# ratio>1 == pappus tip > achene max)
Pap.tip.ach.max.width.ratio_f <- traits_mod_f$Pappus_tip_width/traits_mod_f$Achene_max_width_.eyeballed.
Pap.tip.ach.max.width.ratio_s <- traits_mod_s$Pappus_tip_width/traits_mod_s$Achene_max_width_.eyeballed.

# Finding max and min of the two calculated ratios
Pap.tip.ach.max.width.ratio.max <- pmax(Pap.tip.ach.max.width.ratio_f, 
                                        Pap.tip.ach.max.width.ratio_s)
Pap.tip.ach.max.width.ratio.min <- pmin(Pap.tip.ach.max.width.ratio_f, 
                                        Pap.tip.ach.max.width.ratio_s)

# Add pappus tip to pappus base max ratio traits to aggregated trait data
traits_mod_agg$Pap.tip.base.width.ratio_f <- Pap.tip.base.width.ratio_f
traits_mod_agg$Pap.tip.base.width.ratio_s <- Pap.tip.base.width.ratio_s
traits_mod_agg$Pap.tip.base.width.ratio.max <- Pap.tip.base.width.ratio.max
traits_mod_agg$Pap.tip.base.width.ratio.min <- Pap.tip.base.width.ratio.min

# Add pappus tip to achene max ratio traits to aggregated trait data

traits_mod_agg$Pap.tip.ach.max.width.ratio_f <- Pap.tip.ach.max.width.ratio_f
traits_mod_agg$Pap.tip.ach.max.width.ratio_s <- Pap.tip.ach.max.width.ratio_s
traits_mod_agg$Pap.tip.ach.max.width.ratio.max <- Pap.tip.ach.max.width.ratio.max
traits_mod_agg$Pap.tip.ach.max.width.ratio.min <- Pap.tip.ach.max.width.ratio.min

## -------------
## Length and Length Ratio Traits -------

# Average Achene Length
Ach.len.lin.mean <- apply(cbind(traits_mod_f$Achene_length_lin, traits_mod_s$Achene_length_lin), 1, mean)

# Average Pappus Length
Pap.perp.max.len.mean <- apply(cbind(traits_mod_f$Pappus_perp_max_length, traits_mod_s$Pappus_perp_max_length), 1, mean)

# Add to aggregated trait data
traits_mod_agg$Ach.len.lin.mean <- Ach.len.lin.mean
traits_mod_agg$Pap.perp.max.len.mean <- Pap.perp.max.len.mean

# Creating ratios of pappus perp length to achene length (linear) - made from the 
# average of the front and side measurements
Pap.len.ach.len.ratio <- traits_mod_agg$Pap.perp.max.len.mean/traits_mod_agg$Ach.len.lin.mean

# Add to aggregated trait data
traits_mod_agg$Pap.len.ach.len.ratio <- Pap.len.ach.len.ratio

## -------------
## Angle Traits -------

# Two-Point Angles
traits_mod_f$Achene_twopt_angle_dev <- 180-traits_mod_f$Achene_twopt_angle
traits_mod_s$Achene_twopt_angle_dev <- 180-traits_mod_f$Achene_twopt_angle

# Add in the two point angle in case that matters? This is a proxy for variation 
# in the shape of the cypsela in general
traits_mod_agg$Achene.twopt.angle.dev_f <- traits_mod_f$Achene_twopt_angle_dev
traits_mod_agg$Achene.twopt.angle.dev_s <- traits_mod_s$Achene_twopt_angle_dev

# Max Angle
traits_mod_agg$Achene.twopt.angle.dev.max <- 
  pmax(traits_mod_agg$Achene.twopt.angle.dev_f, traits_mod_agg$Achene.twopt.angle.dev_s)

traits_mod_agg$Achene.twopt.angle.dev.min <- 
  pmin(traits_mod_agg$Achene.twopt.angle.dev_f, traits_mod_agg$Achene.twopt.angle.dev_s)


## -------------
## Area and Area Ratio Traits ------

# Defining pi for area calculations
pi = 3.14159265359

# Area Traits
pappus_tip_area <- pi*(traits_mod_f$Pappus_tip_width/2)*(traits_mod_s$Pappus_tip_width/2)
pappus_base_area <- pi*(traits_mod_f$Pappus_base_width/2)*(traits_mod_s$Pappus_base_width/2)
cypsela_max_area <- pi*(traits_mod_f$Achene_max_width_.eyeballed./2)*(traits_mod_s$Achene_max_width_.eyeballed./2)

# Area Ratio Traits
pap.tip.cyp.max.area.ratio <- pappus_tip_area/cypsela_max_area
pap.tip.pap.base.area.ratio <- pappus_tip_area/pappus_base_area
traits_mod_agg$pappus_tip_area <- pappus_tip_area
traits_mod_agg$cypsela_max_area <- cypsela_max_area
traits_mod_agg$pap.tip.cyp.max.area.ratio <- pap.tip.cyp.max.area.ratio
traits_mod_agg$pap.tip.pap.base.area.ratio <- pap.tip.pap.base.area.ratio


## -------------
## Eccentricity and Surface Area Traits ------

# Pappus Tip Eccentricity
 
# Find the semi-major and minor axes for each pappus tip ellipse
major_ax_tip <- apply(cbind(traits_mod_f$Pappus_tip_width, traits_mod_s$Pappus_tip_width), 1, max)/2
minor_ax_tip <- apply(cbind(traits_mod_f$Pappus_tip_width, traits_mod_s$Pappus_tip_width), 1, min)/2

# Calculate eccentricity
pap_tip_ecc <- sqrt(1-((minor_ax_tip^2)/(major_ax_tip^2)))

# Pappus Base Eccentricity
# Find the semi-major and minor axes for each pappus base ellipse (this is the 
# "radius" in a given direction) 
major_ax_base <- apply(cbind(traits_mod_f$Pappus_base_width, traits_mod_s$Pappus_base_width), 1, max)/2
minor_ax_base <- apply(cbind(traits_mod_f$Pappus_base_width, traits_mod_s$Pappus_base_width), 1, min)/2

# Calculate eccentricity
pap_base_ecc <- sqrt(1-((minor_ax_base^2)/(major_ax_base^2)))

# Add to aggregated traits df
traits_mod_agg$pap.tip.ecc <- pap_tip_ecc
traits_mod_agg$pap.base.ecc <- pap_base_ecc


## -------------
## Surface Area Traits ------

# Paraboloid Surface Area
ach.cone.data <- data.frame("a" = major_ax_base, "b" = minor_ax_base, "h" = traits_mod_agg$Ach.len.lin.mean)

# Calculating the surface area of an elliptic parabaloid in R 
# found here: 
# https://stackoverflow.com/questions/22719335/how-to-implement-surface-calculation-of-a-elliptical-cone-and-elliptical-parabol

# Calculate the (lateral) surface area of an elliptic paraboloid 

parabA <- function(a,b,h) {
  require(cubature)
  integrand <- function(x,a,b){
    u <- x[1]
    v <- x[2]
    E <- 1+((a*cos(v))^2 + (b*sin(v))^2)/(4*u)
    F1 <- (b^2 - a^2)*sin(2*v)/4
    G <- u*((a*sin(v))^2+(b*cos(v))^2)
    return(sqrt(E*G-F1^2))
  }
  adaptIntegrate(integrand, lowerLimit=c(0,0), upperLimit=c(h,2*pi), a=a, b=b)$integral
}

# Make the data vector 
parabSA <- mapply(parabA, a = ach.cone.data$a, b = ach.cone.data$b, h = ach.cone.data$h)

# Cone Surface Area
# use the defined elliptic cone formulae here: http://mathworld.wolfram.com/EllipticCone.html 
# and then modified the code above.

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

# Make the data vector
coneSA <- mapply(coneA, a = ach.cone.data$a, b = ach.cone.data$b, h = ach.cone.data$h)

# Add traits to aggregated data
traits_mod_agg$parabSA <- parabSA
traits_mod_agg$coneSA <- coneSA

# End product preview
traits_mod_agg[1:3,]


### ------------------
### Split by Species -------

# Split trait data frame by species
traits_mod_agg.cal <- arrange(filter(traits_mod_agg, Species == "L. californica"), Tube.ID)
traits_mod_agg.fre <- arrange(filter(traits_mod_agg, Species == "L. fremontii"), Tube.ID)
traits_mod_agg.gla <- arrange(filter(traits_mod_agg, Species == "L. glaberrima"), Tube.ID)


traits_mod_agg.cal.mod <- traits_mod_agg.cal
traits_mod_agg.fre.mod <- traits_mod_agg.fre
traits_mod_agg.gla.mod <- traits_mod_agg.gla


# Split distance data frame by species
dists2_main.cal <- arrange(filter(dists2_main, Species == "californica"), Tube.ID)
dists2_main.fre <- arrange(filter(dists2_main, Species == "fremontii"), Tube.ID)
dists2_main.gla <- arrange(filter(dists2_main, Species == "glaberrima"), Tube.ID)
