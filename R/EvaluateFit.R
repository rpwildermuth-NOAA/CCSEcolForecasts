###########################################################################################################
# Evaluate the fit of an SDM on withheld test data
# Steps are to 1) build an SDM using a specifified number of training years, and
# 2) evaluate SDM skill on the year/s after the training data ends (i.e. X year forecast)
# Step 2 includes assessing skill across multiple spatial and temporal resolutions
# Now has ability to also look at persistence forecast (i.e. environmental predictors from previous year)
# For now is set up to work with a presence/absence (binomial) SDM
# Contact bmuhling@ucsc.edu
###########################################################################################################

library(mgcv)
library(here)
library(tidyr)
library(lubridate)
library(pROC)
library(ggplot2)
library(scales)
library(dplyr)
library(gbm)
library(dismo)
library(viridis)

###########################################################################################################
# Load the helper functions you will need
source("./R/scoreSDM.R")
source("./R/testSkillSDM.R")
source("./R/runDiagnostics.R")

###########################################################################################################
# Load biological observation data with environmental covariates extracted
# This is now updated to include variables extracted at a lag of 1 year, which are used to calculate
# persistence skill (variables labeled with "_lag1")
# Is assumed that the columns "lon", "lat", and "date" are present and complete
obs <- readRDS("./data/combinedBioObs_envExtracted.rds")
# If observations don't have a "year", "month" and "quarter" columns, add them. Will need a "date" column to get them
if(!"year" %in% colnames(obs)) {
  obs$year <- year(obs$date)  
}
if(!"month" %in% colnames(obs)) {
  obs$month <- month(obs$date)  
}
if(!"quarter" %in% colnames(obs)) {
  obs$quarter <- quarter(obs$date)  
}

# I'm missing predictor variables for 2024, so I'm trimming that year off
obs <- subset(obs, year <= 2023)
# Define the target variable: here presence/absence of anchovy
obs$pa <- obs$anchPA

###########################################################################################################
# Define some parameters for building the SDM/s
# Note that for now the user supplies the SDM-specific parameters. Optimizing parameters for an SDM (esp. a BRT) 
# when you don't have good information on dataset size, number of variables, prevalence rate etc. is difficult
# This could be improved in future iterations of this code!
sdmType <- "brt" # Can currently be "gam" or "brt"
k <- 4 # Number of knots for a GAM. Will be ignored if not building a GAM
tc <- 3 # tree complexity for a BRT. Will be ignored if not building a BRT
lr <- 0.02 # learning rate for a BRT. Will be ignored if not building a BRT
max.trees <- 10000 # max trees for a BRT. Will be ignored if not building a BRT
varNames <- c("sst", "ild", "sst_sd", "ssh_sd", "logChl", "logEKE", "distLand", "anchssb") # Names of predictors in the SDM
targetName <- "pa" # Target variable for the SDM
aucCutoff <- 10 # If less observations than this cutoff within a month/season etc., AUC will not be calculated

###########################################################################################################
# Function to 1) build an SDM and predict to X years of withheld data using scoreSDM, and 
# 2) calculate SDM skill at various spatial/temporal aggregation using a specified subset of the 
# observational dataset for training. "noYrs" is the number of years to use for SDM training
peelSDM <- function(obs, sdmType, varNames, targetName, k, tc, lr, max.trees, noYrs, 
                    yrsToForecast, includePersistence) {
  # Define subObs: the subset of the complete data to be used for training, plus forecast/testing year/s
  yrs <- unique(sort(obs$year))
  # Define years to include in subObs. At present, training data always start at the earliest year available,
  # this could be adjusted
  yrsToInclude <- yrs[1: (noYrs + yrsToForecast)] 
  
  # Catch if you try to forecast beyond the available data
  if(length(yrs) < length(yrsToInclude)) {
    stop(print("Attempting to forecast beyond available data"))
  }
  
  subObs <- subset(obs, year %in% yrsToInclude)
  
  # Call scoreSDM to build an SDM using subObs as training and test data
  mod <- scoreSDM(subObs = subObs, sdmType = sdmType, varNames = varNames,
                       targetName = targetName, k = k, tc = tc, lr = lr, max.trees = max.trees, 
                       yrsToForecast = yrsToForecast, includePersistence = includePersistence)
  # Call testSkillSDM to assess the skill of the X-year forecast for this SDM
  sdmSkill <- testSkillSDM(mod = mod, targetName = targetName, aucCutoff = aucCutoff, usePersistence = FALSE)
  # Add noYrs and sdmType to output
  sdmSkill$noYrsTrain <- noYrs
  sdmSkill$terminalYr <- max(subObs$year) - yrsToForecast 
  sdmSkill$sdmType <- sdmType
  
  # If you want to assess persistence skill, add a new list element showing sdmSkill for those predictions
  if(includePersistence == TRUE) {
    sdmSkillPers <- testSkillSDM(mod = mod, targetName = targetName, aucCutoff = aucCutoff, usePersistence = TRUE)
    # Add noYrs and sdmType to output
    sdmSkillPers$noYrsTrain <- noYrs
    sdmSkillPers$terminalYr <- max(subObs$year) - yrsToForecast 
    sdmSkillPers$sdmType <- sdmType
  }
  
  # Return some simple diagnostics
  evalOutputs <- runDiagnostics(mod = mod, max.trees = max.trees)
  
  # Return results: at present, returning as a list
  if(includePersistence == TRUE) {
    out <- list("sdm" = mod$sdm, "test" = mod$test, "modelEval" = evalOutputs, "sdmSkill" = sdmSkill, "sdmSkillPers" = sdmSkillPers)
  } else {
    out <- list("sdm" = mod$sdm, "test" = mod$test, "modelEval" = evalOutputs, "sdmSkill" = sdmSkill)
  }
  return(out)
}

###########################################################################################################
# An example: 5-year forecast skill for anchovy SDM where training data include between 10 and 15 years of data 
yrsToTrain <- seq(10, 15) 
yrsToForecast <- 5
# Optional: do we want to assess persistence skill? That is: can last year's environmental conditions
# predict this years species distributions?
includePersistence <- TRUE
output <- vector(mode = "list", length = length(yrsToTrain)) # We will save results to a list (a list of lists)
for (j in 1:length(output)) {
  output[[j]] <- peelSDM(obs = obs, sdmType = sdmType, varNames = varNames, targetName = targetName,
                         k = k, tc = tc, lr = lr, max.trees = max.trees, noYrs = yrsToTrain[j], 
                         yrsToForecast = yrsToForecast, includePersistence = includePersistence)
}

###########################################################################################################
# Some quick plots of SDM skill with varying lengths of training data, and varying levels of spatial/temporal aggregation
# for each forecasting horizon separately
# First reshape the output list so is easier to work with
suppressWarnings(rm(sdmSkill1, sdmSkillAll))
for(b in 1:length(output)) { 
  sdmSkill1 <- output[[b]]$sdmSkill
  if(exists("sdmSkillAll")) {
    sdmSkillAll <- rbind(sdmSkillAll, sdmSkill1)
  } else {
    sdmSkillAll <- sdmSkill1
  }
}

# Aggregate so skills are averaged across different months/seasons. "year" is the test year
sdmSkillAllAgg <- aggregate(auc ~ time + space + noYrsTrain + terminalYr + year, sdmSkillAll, FUN = mean, na.rm = TRUE)
# Add a field converting year to forecast time horizon
sdmSkillAllAgg$forecastHorizon <- sdmSkillAllAgg$year - sdmSkillAllAgg$terminalYr

# Plot showing mean AUC skill across all terminal years by spatial and temporal aggregation and forecast horizon
ggplot(sdmSkillAllAgg) + 
  geom_boxplot(aes(x = forecastHorizon, y = auc, group = factor(forecastHorizon), fill = forecastHorizon)) +
  scale_fill_viridis("Forecast \nHorizon", option = "mako") + xlab("Forecast Horizon (Years)") + ylab("Forecast AUC") + 
  scale_y_continuous(limits = c(0, 1), oob = rescale_none) + theme_bw() + facet_grid(time ~ space, scales = "free")

###########################################################################################################
# An additional plot showing real/contemporary skill vs persistence skill, if you also calculated and returned those results
# Reshape the persistence skill outputs
suppressWarnings(rm(sdmSkillPers1, sdmSkillPersAll))
for(b in 1:length(output)) { 
  sdmSkillPers1 <- output[[b]]$sdmSkillPers
  if(exists("sdmSkillPersAll")) {
    sdmSkillPersAll <- rbind(sdmSkillPersAll, sdmSkillPers1)
  } else {
    sdmSkillPersAll <- sdmSkillPers1
  }
}

# Aggregate so skills are averaged across different months/seasons. "year" is the test year
sdmSkillPersAllAgg <- aggregate(auc ~ time + space + noYrsTrain + terminalYr + year, sdmSkillPersAll, FUN = mean, na.rm = TRUE)
# Add a field converting year to forecast time horizon
sdmSkillPersAllAgg$forecastHorizon <- sdmSkillPersAllAgg$year - sdmSkillPersAllAgg$terminalYr
# Combine with skill assessment from contemporary data
sdmSkillAllAgg$forecastType <- "contemporary"
sdmSkillPersAllAgg$forecastType <- "persistence"
sdmSkillBoth <- rbind(sdmSkillAllAgg, sdmSkillPersAllAgg)
sdmSkillBoth$horizonType <- interaction(sdmSkillBoth$forecastHorizon, sdmSkillBoth$forecastType)

# A plot showing skill by aggregation level and forecast horizon, now comparing contemporary prediction vs. persistence
ggplot() + 
  geom_boxplot(data = sdmSkillBoth, aes(x = forecastHorizon, y = auc, group = factor(horizonType), fill = forecastHorizon)) +
  scale_fill_viridis("Forecast \nHorizon", option = "mako") + xlab("Forecast Horizon (Years)") + ylab("Forecast AUC") + 
  scale_y_continuous(limits = c(0, 1), oob = rescale_none) + theme_bw() + facet_grid(time ~ space)

###########################################################################################################
# Save output. Build a more informative filename depending on your work!
saveRDS(output, file = paste0("./outputs/", sdmType, "_anchovy_", targetName, "forecast", yrsToForecast, "years.rds")) 
