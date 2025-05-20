###########################################################################################################
# Extract environmental predictors at locations of biological sampling/observation points
# Now operates as a function which intakes an observational dataset, a list of variables to extract
# and a lag. Lag = 1 will extract variables from 1 year ago: used to calculate persistence forecast skill
# Lag = 0 extracts environmental variables at the date that sampling occurred
# Contact bmuhling@ucsc.edu
###########################################################################################################

library(here)
library(sf)
library(lubridate)
library(dplyr)
library(tidyr)

###########################################################################################################
# Load helper functions to extract environmental predictors
# These rely on environmental netcdfs, which are stored in different places!
source("./R/getROMS.R") # ROMS ocean model variables
# source("./R/getMOM6_gridded.R") # MOM6 ocean model variables
source("./R/getCMEMS_l4chl.R") # CMEMS daily L4 chl
source("./R/getDistLand.R") # Distance to nearst land based on coast shp
source("./R/getBathym.R") # ETOPO bathymetry
# Load datasets: these currently live in project ./data
coast <- sf::read_sf("./data/EPOCoast60_noGI.shp") # The coast shp
anch <- read.csv("./data/Hinchliffe_CSNA_timeseries_19652023.csv") # Anchovy SSB from Hincliffe et al. 2025

# Load biological observations
# They should have columns labelled "lon" (longitude), "lat" (latitude) and "date" (date in R format)
# (If you load a csv, you may need to re-format the date)
obs <- readRDS("./data/combinedBioObs.rds")

# Different SDMs use many different variables
# Rather than attempting to pass in a list of variables that covers all possibilities, 
# this script extracts a specific set of variables for an anchovy SDM using a series of external functions 
# This could be adjusted for another species/set of predictors
# Note that "desired.diameter" defines a box around an observation within which to extract environmental variable
# e.g. desired.diameter = 0.7 returns values within a 0.7x0.7 degree box
# Here I will extract these variables for an anchovy SDM:
# ROMs: ild, sst, logEKE (really TKE), sst_sd (0.7 degrees), ssh_sd (0.7 degrees)
# MOM6 (gridded): sst, surface chl 
# CMEMS: surface chl
# SSB (Hinchliffe et al.)
# Distance from nearest land (relative to a coastline shapefile)
###########################################################################################################

# To get a "persistence" forecast, we extract environmental variables from the year before the sampling date
# lag = 0 in the function extracts variables at the sampling date, lag = 1 extracts from the year before
# This is clunky but here I'm running the function twice to get lag = 0 and lag = 1, then joining the outputs
# Careful with lag > 0, could potentially give dates outside the temporal range of environmental datasets
envExtractLag0 <- extractEnvVars(obs = obs, lag = 0)
envExtractLag1 <- extractEnvVars(obs = obs, lag = 1)
# Adjust colnames: clunky
colnames(envExtractLag1) <- paste0(colnames(envExtractLag1), "_lag1") # Adjust with lag number
# Figure out col index where environmental covariates start
startCol <- ncol(obs) + 1
envExtract <- cbind(envExtractLag0, envExtractLag1[, startCol: ncol(envExtractLag1)])
# Save the output
saveRDS(envExtract, "./data/combinedBioObs_envExtracted.rds")

###########################################################################################################
# The extraction function
extractEnvVars <- function(obs, lag) {
  # Adjust the date depending on the lag
  if (lag == 0) {
    obs$date <- obs$date
  } else if (lag > 0) {
    obs$date <- obs$date - years(x = lag)
  }
  
  # Extract the environmental variables
  # ROMS variables
  print("Extracting ROMS Variables")
  ildMean <- getROMS(points = obs, varName = "ild", desired.diameter = 0.1, func = "mean", 
                     histPath = "F:/roms/hist", nrtPath = "F:/roms/nrtComplete")
  sstMean <- getROMS(points = obs, varName = "sst", desired.diameter = 0.1, func = "mean", 
                     histPath = "F:/roms/hist", nrtPath = "F:/roms/nrtComplete")
  sstSD <- getROMS(points = obs, varName = "sst", desired.diameter = 0.7, func = "sd", 
                   histPath = "F:/roms/hist", nrtPath = "F:/roms/nrtComplete")
  sshSD <- getROMS(points = obs, varName = "ssh", desired.diameter = 0.7, func = "sd", 
                   histPath = "F:/roms/hist", nrtPath = "F:/roms/nrtComplete")
  # Current vectors for EKE
  suMean <- getROMS(points = obs, varName = "su", desired.diameter = 0.1, func = "mean", 
                    histPath = "F:/roms/hist", nrtPath = "F:/roms/nrtComplete")
  svMean <- getROMS(points = obs, varName = "sv", desired.diameter = 0.1, func = "mean", 
                    histPath = "F:/roms/hist", nrtPath = "F:/roms/nrtComplete")
  # Join to observations
  obs$ild <- ildMean$ild_mean_0.1
  obs$sst <- sstMean$sst_mean_0.1
  obs$sst_sd <- sstSD$sst_sd_0.7
  obs$ssh_sd <- sshSD$ssh_sd_0.7
  obs$eke <- ((suMean$su_mean_0.1 ^ 2) + (svMean$sv_mean_0.1 ^ 2)) / 2
  obs$logEKE <- log(obs$eke)
  
  ###########################################################################################################
  # MOM6 variables
  # print("Extracting MOM6 Variables")
  ### Note daily surface chl not currently available from CEFI portal for regular grid ###
  ### Will update here when it is available ###
  
  ###########################################################################################################
  # CMEMS level 4 surface chlorophyll
  print("Extracting CMEMS Surface Chl")
  chlMean <- getCMEMS_l4chl(points = obs, desired.diameter = 0.1, func = "mean", nc.path = "F:/")
  # Join to observations
  obs$chl <- chlMean$chl_mean_0.1
  obs$logChl <- log(obs$chl)
  
  ###########################################################################################################
  # Distance from shore: requires a coast shapefile. Can be slow
  print("Extracting Distance From Land")
  distLand <- getDistLand(points = obs, coast = coast)
  # Join to observations
  obs$distLand <- distLand$distLand
  
  ###########################################################################################################
  # Bathymetry
  print("Extracting Water Depth")
  bathym <- getBathym(points = obs, desired.diameter = 0.1, func = "mean")
  obs$bathym <- bathym$bathym_mean_0.1
  
  ###########################################################################################################
  # Anchovy SSB (from Hinchliffe et al. 2025)
  anch$year <- anch$Model.Y.S
  anch$year <- as.numeric(gsub("-1|-2", "", anch$year))
  anch <- subset(anch, !is.na(anch$SSB..mt.))
  anch$anchssb <- anch$SSB..mt.
  # Join to observations
  obs$year <- year(obs$date)
  obs <- left_join(obs, anch[c("year", "anchssb")], by = "year")
  
  ###########################################################################################################
  # Save observations with environmental variables extracted
  # saveRDS(obs, "./data/combinedBioObs_envExtracted.rds")
  # Return the DF
  return(obs)
}
