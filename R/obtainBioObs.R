#########################################################################################################
# How to build an anchovy SDM: step 1:
# Collate biological observations
# Two sources of data: one is the NOAA SWFSC CPS trawl cruises. Data are freely available on ERDDAP
# The other is the Columbia River Predator surveys. These data are available by request only from the
# NWFSC, and so are loaded from a local folder
#########################################################################################################

library(rerddap)
library(tidyr)
library(lubridate)
library(readxl)
library(here)

dataPath <- here::here("cefi", "reforecast", "data")

#########################################################################################################
######################### Observations from CPS Trawl survey ############################################
#########################################################################################################
# Download the latest SWFSC CPS trawl observations from ERDDAP
dataInfo <- info('FRDCPSTrawlLHHaulCatch', url = 'https://coastwatch.pfeg.noaa.gov/erddap/')
# Show the dataset columns (also see https://coastwatch.pfeg.noaa.gov/erddap/tabledap/FRDCPSTrawlLHHaulCatch.html)
cols <- info('FRDCPSTrawlLHHaulCatch')$variables
# Download the data
cps <- tabledap(dataInfo, fields = cols$variable_name)

# Process fields 
# NA values for subsamples mean zeroes
cps$subsample_count <- ifelse(is.na(cps$subsample_count), 0, cps$subsample_count)
cps$subsample_weight <- ifelse(is.na(cps$subsample_weight), 0, cps$subsample_weight)
# Use numbers and weights to show presence/absence
cps$pres <- ifelse(cps$subsample_count > 0 | cps$subsample_weight > 0, 1, 0)
# Get mean longitude and latitude
cps$lon <- rowMeans(cps[c("longitude", "stop_longitude")])
cps$lat <- rowMeans(cps[c("latitude", "stop_latitude")])
# plot(cps$lon, cps$lat) # Quick for any outlier locations
# Add a date column
cps$date <- as.Date(cps$time)
# Remove a few daytime samples
cps$hr <- hour(cps$time - hours(7)) # times originally in UTC
cps$dn <- ifelse(cps$hr < 6.1 | cps$hr > 17.9, "night", "day")
cps <- subset(cps, dn == "night")
# Convert to wide format
cpsMatrix <- pivot_wider(cps, names_from = scientific_name, values_from = pres, 
                         id_cols = c(cruise, haul, lon, lat, date, dn), values_fill = list(values = 0))
# Convert NA to 0 
cpsMatrix[is.na(cpsMatrix)] <- 0

# Add a few other useful columns
cpsMatrix$survey <- "cps"

# Subset to just anchovy
cpsAnch <- cpsMatrix[c("survey", "cruise", "haul", "lon", "lat", "date", "Engraulis mordax")]
colnames(cpsAnch)[ncol(cpsAnch)] <- "anchPA"

#########################################################################################################
######################### Observations from Columbia River Predator Surveys##############################
#########################################################################################################
predator <- read_excel(here::here(dataPath, "predator", "Muhling_PredatorStudy_Forage_CPUE_7.30.20.xlsx"))
# Format columns as needed
predator$lon <- predator$`Start Long (decimal degrees)`
predator$lat <- predator$`Start Lat (decimal degrees)`
predator$date <- as.Date(predator$`Start Time`)
predator$hr <- hour(predator$`Start Time`)
predator$cruise <- year(predator$date) # No cruise identifier, and haul is unique, so just using year
predator$haul <- predator$`Haul #`
predator$dn <- ifelse(predator$hr < 6.1 | predator$hr > 17.9, "night", "day")
# Set presences/absences
predator$anchPA <- ifelse(predator$`Northern anchovy` > 0, 1, 0)
# Remove day samples, remove "non-valid" samples, remove #4 rope trawl samples (per Cheryl Mogran, NWFSC)
predator <- subset(predator, predator$dn == "night" & predator$valid == 1 & predator$`net type` != "#4 rope trawl")
predator$survey <- "predator"
predator <- predator[c("survey", "cruise", "haul", "lon", "lat", "date", "anchPA")]

#########################################################################################################
######################### Join the two datasets and save them ###########################################
#########################################################################################################
allData <- rbind(cpsAnch, predator)
# Add a unique ID
allData$id <- seq(1:nrow(allData))
saveRDS(allData, "./data/combinedBioObs.rds")
