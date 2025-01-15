###########################################################################################################
# Evaluate the fit of an SDM on withheld test data
# Testing withholding the last 1 - 10 years of available time-series
# This example shows a GAM for anchovy presence/absence
# Contact Barbara.Muhling@noaa.gov
###########################################################################################################

library(mgcv)
library(here)
library(tidyr)
library(lubridate)
library(pROC)
library(ggplot2)
library(scales)

###########################################################################################################
# Load biological observation data with environmental covariates extracted
# (Mine live outside the project)
parent <- here() %>% dirname()
dataPath <- here::here(parent, "cefi", "reforecast", "data")
obs <- readRDS(here::here(dataPath, "combinedBioObs_envExtracted.rds"))
# If observations don't have a "year" column, add one
# obs$year <- year(obs$date)
# I'm missing predictor ariables for 2024, so trimming that year off
obs <- subset(obs, year <= 2023)
# Define the target variable
obs$pa <- obs$anchPA

###########################################################################################################
# Make a small function to build a GAM withholding the last X years of data
# Can set maximum k if desired
buildGAM <- function(obs, withhold, k) {
  allYrs <- aggregate(pa ~ year, obs, FUN = length)
  # Define the last X years (will not include years with no observations, e.g. 2020)
  testYrs <- allYrs$year[(nrow(allYrs) - (withhold - 1)) : nrow(allYrs)]
  test <- subset(obs, year %in% testYrs)
  train <- subset(obs, !year %in% testYrs)
  
  # Now build a GAM. I am keeping SSB partial response linear
  gam1 <- gam(pa ~ s(sst, k = k) + s(ild, k = k) + s(sst_sd, k = k) + s(ssh_sd, k = k) + s(logChl, k = k) + 
                s(logEKE, k = k) + s(distLand, k = k) + s(anchssb, k = 3),
              data = train, family = "binomial", method = "REML", select = TRUE)
  # summary(gam1)
  # plot(gam1, scale = 0, pages = 1)
  
  # Calculate the AUC for training and testing datasets
  aucs <- data.frame("trainAUC" = NA, "testAUC" = NA)
  train$pred <- predict(gam1, train, type = "response")
  test$pred <- predict(gam1, test, type = "response")
  aucs$trainAUC <- auc(train$pa, train$pred, direction = "<", quiet = TRUE)
  aucs$testAUC <- auc(test$pa, test$pred, direction = "<", quiet = TRUE) ###
  
  # Return everything
  # out <- list("gam1" = gam1, "train" = train, "test" = test, "AUCs" = aucs)
  # Or just return AUCs for now
  out <- aucs
  return(out)
}

###########################################################################################################
# Loop through to get train and test AUCs from withholding last 1 - 10 years of observations
allAUCs <- data.frame(matrix(nrow = 20, ncol = 3))
colnames(allAUCs) <- c("yearsWithheld", "trainAUC", "testAUC")
for (i in 1:nrow(allAUCs)) {
  allAUCs[i, 1] <- i
  allAUCs[i, 2:3] <- buildGAM(obs = obs, withhold = i, k = 4)
}

# Quick plot
aucsLong <- pivot_longer(allAUCs, cols = c("trainAUC", "testAUC"), names_to = "testTrain", values_to = "auc")
ggplot(aucsLong) + geom_bar(aes(x = yearsWithheld, y = auc, fill = testTrain), stat = "identity", position = "dodge") + 
  scale_y_continuous(limits = c(0.5, 0.9), oob = rescale_none) + theme_bw()
