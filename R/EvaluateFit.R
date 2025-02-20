###########################################################################################################
# Evaluate the fit of an SDM on withheld test data
# Steps are to 1) build an SDM using a specifified number of training years, and
# 2) evaluate SDM skill on the year after the training data ends (i.e. 1 year forecast)
# Step 2 includes assessing skill across multiple spatial and temporal resolutions
# Contact Barbara.Muhling@noaa.gov
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

###########################################################################################################
# Load the helper functions you will need
source("./R/scoreSDM.R")
source("./R/testSkillSDM.R")

###########################################################################################################
# Load biological observation data with environmental covariates extracted
# Is assumed that the columns "lon", "lat", and "date" are present and complete
# parent <- here() %>% dirname()
# dataPath <- here::here(parent, "cefi", "reforecast", "data")
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
sdmType <- "gam" # Can currently be "gam" or "brt"
k <- 4 # Number of knots for a GAM. Will be ignored if not building a GAM
tc <- 3 # tree complexity for a BRT. Will be ignored if not building a BRT
lr <- 0.025 # learning rate  for a BRT. Will be ignored if not building a BRT
varNames <- c("sst", "ild", "sst_sd", "ssh_sd", "logChl", "logEKE", "distLand", "anchssb") # Names of predictors in the SDM
targetName <- "pa" # Target variable for the SDM
aucCutoff <- 10 # If less obervations than this cutoff within a month/season etc., AUC will not be calculated

###########################################################################################################
# Function to 1) build an SDM and predict to 1 year of withheld data using scoreSDM, and 
# 2) calculate SDM skill at various spatial/temporal aggregation using a specified subset of the 
# observational dataset for training. "noYrs" is the number of years to use for SDM training
peelSDM <- function(obs, sdmType, varNames, targetName, k, tc, lr, noYrs) {
  # Define subObs: the subset of the complete data to be used for training, plus 1 forecast/testing year
  yrs <- unique(sort(obs$year))
  # Define years to include in subObs. At present, training data always start at the earliest year available,
  # this could easily be adjusted
  yrsToInclude <- yrs[1: (noYrs + 1)] 
  subObs <- subset(obs, year %in% yrsToInclude)
  
  # Call scoreSDM to build an SDM using subObs as training and test data
  mod <- scoreSDM(subObs = subObs, sdmType = sdmType, varNames = varNames,
                                               targetName = targetName, k = k, tc = tc, lr = lr)
  # Call testSkillSDM to assess the skill of the 1-year forecast for this SDM
  sdmSkill <- testSkillSDM(mod = mod, targetName = targetName, aucCutoff = aucCutoff)
  # Add noYrs and sdmType to output
  sdmSkill$noYrsTrain <- noYrs
  sdmSkill$terminalYr <- max(subObs$year) - 1
  sdmSkill$sdmType <- sdmType
  
  # Return results: at present, returning mod and sdmSkill as a list
  out <- list("sdm" = mod$sdm, "test" = mod$test, "sdmSkill" = sdmSkill)
  return(out)
}

###########################################################################################################
# An example: 1-year forecast skill for anchovy SDM where training data include between 10 and 20 years of data 
yrsToTrain <- seq(10, 24) 
output <- vector(mode = "list", length = length(yrsToTrain)) # We will save results to a list (a list of lists)
for (j in 1:length(output)) {
  output[[j]] <- peelSDM(obs = obs, sdmType = sdmType, varNames = varNames, targetName = targetName,
                         k = k, tc = tc, lr = lr, noYrs = yrsToTrain[j])
}

###########################################################################################################
# A quick plot of SDM skill with varying lengths of training data, and varying levels of spatial/temporal aggregation
suppressWarnings(rm(sdmSkill1, sdmSkillAll))
for(k in 1:length(output)) { 
  sdmSkill1 <- output[[k]]$sdmSkill
  if(exists("sdmSkillAll")) {
    sdmSkillAll <- rbind(sdmSkillAll, sdmSkill1)
  } else {
    sdmSkillAll <- sdmSkill1
  }
}

# Aggregate so skills are averaged across different months/seasons
sdmSkillAllAgg <- aggregate(auc ~ time + space + noYrsTrain + terminalYr, sdmSkillAll, FUN = mean, na.rm = TRUE)

# Plot
ggplot(sdmSkillAllAgg) + geom_bar(aes(x = terminalYr, y = auc, group = time, fill = time), 
                                  stat = "identity", position = position_dodge(preserve = "single")) +
    scale_y_continuous(limits = c(0.5, 1), oob = rescale_none) + theme_bw() + facet_wrap(~ space, ncol = 1)

# Save output. Build a more informative filename depending on your work!
saveRDS(output, file = paste0("./outputs/", sdmType, "_anchovy_", targetName, ".rds")) 