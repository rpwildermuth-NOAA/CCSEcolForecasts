###########################################################################################################
# Extract environmental predictors at locations of biological sampling/observation points
# Contact Barbara.Muhling@noaa.gov
###########################################################################################################

library(here)
library(sf)
library(lubridate)
library(dplyr)
library(tidyr)

# # Define the path where biological observations are (now contained inside the project)
# parent <- here() %>% dirname()
# dataPath <- here::here(parent, "cefi", "reforecast", "data")

###########################################################################################################
# Load biological observations
# They should have columns labelled "lon" (longitude), "lat" (latitude) and "date" (date in R format)
#(If you load a csv, you may need to re-format the date)
obs <- readRDS("./data/combinedBioObs.rds")

###########################################################################################################
# Use a series of external functions to extract desired environmental predictors
# Note that "desired.diameter" defines a box around an observation within which to extract environmental variable
# e.g. desired.diameter = 0.7 returns values within a 0.7x0.7 degree box
# Here I will extract variables for an anchovy SDM:
# ROMs: ild, sst, logEKE (really TKE), sst_sd (0.7 degrees), ssh_sd (0.7 degrees)
# MOM6 (gridded): sst, surface chl 
# CMEMS: surface chl
# SSB (Hinchliffe et al.)
# Distance from nearest land (relative to a coastline shapefile)
###########################################################################################################
# ROMS variables
source("./R/getROMS.R")
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
source("./R/getMOM6.R")
### Note daily surface chl not currently available from CEFI portal for regular grid ###
### Will update here when it is available ###

###########################################################################################################
# CMEMS level 4 surface chlorophyll
source("./R/getCMEMS_l4chl.R")
chlMean <- getCMEMS_l4chl(points = obs, desired.diameter = 0.1, func = "mean", nc.path = "F:/")
# Join to observations
obs$chl <- chlMean$chl_mean_0.1
obs$logChl <- log(obs$chl)

###########################################################################################################
# Distance from shore: requires a coast shapefile. Can be slow
coast <- sf::read_sf(here::here(dataPath, "EPOCoast60_noGI.shp"))
source("./R/getDistLand.R")
distLand <- getDistLand(points = obs, coast = coast)
# Join to observations
obs$distLand <- distLand$distLand

###########################################################################################################
# Bathymetry
source("./R/getBathym.R")
bathym <- getBathym(points = obs, desired.diameter = 0.1, func = "mean")
obs$bathym <- bathym$bathym_mean_0.1

###########################################################################################################
# Anchovy SSB (from Hinchliffe et al. 2025)
anch <- read.csv(here::here(dataPath, "Hinchliffe_CSNA_timeseries_19652023.csv"))
anch$year <- anch$Model.Y.S
anch$year <- as.numeric(gsub("-1|-2", "", anch$year))
anch <- subset(anch, !is.na(anch$SSB..mt.))
anch$anchssb <- anch$SSB..mt.
# Join to observations
obs$year <- year(obs$date)
obs <- left_join(obs, anch[c("year", "anchssb")], by = "year")

###########################################################################################################
# Save observations with environmental variables extracted
saveRDS(obs, "./data/combinedBioObs_envExtracted.rds")
