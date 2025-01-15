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
library(dplyr)

###########################################################################################################
# Load biological observation data with environmental covariates extracted
# (Mine live outside the project)
parent <- here() %>% dirname()
dataPath <- here::here(parent, "cefi", "reforecast", "data")
obs <- readRDS(here::here(dataPath, "combinedBioObs_envExtracted.rds"))
# If observations don't have a "year" column, add one
# obs$year <- year(obs$date)
# I'm missing predictor variables for 2024, so trimming that year off
obs <- subset(obs, year <= 2023)
# Define the target variable
obs$pa <- obs$anchPA

###########################################################################################################
# Make a small function to build a GAM, and forecast to the next X (5?) years
buildGAM <- function(obs, terminalYr, k) {
  # Define the training and test forecast year
  train <- subset(obs, year <= terminalYr)
  test <- subset(obs, year > terminalYr & year <= (terminalYr + 5)) # May be less than 5 years if missing years, e.g. 2020
  
  # Now build a GAM. SSB could be further constrained if necessary
  gam1 <- gam(pa ~ s(sst, k = k) + s(ild, k = k) + s(sst_sd, k = k) + s(ssh_sd, k = k) + s(logChl, k = k) + 
                s(logEKE, k = k) + s(distLand, k = k) + s(anchssb, k = k),
              data = train, family = "binomial", method = "REML", select = TRUE)
  # summary(gam1)
  # plot(gam1, scale = 0, pages = 1)
  
  # Calculate the AUC for training and forecast test datasets
  # First calculate how many forecast years we have (could be some years missing)
  testYrs <- sort(unique(test$year))
  # Now construct a dataframe to receive AUC values
  forecastLabels <- testYrs - terminalYr
  forecastLabels <- paste0(forecastLabels, "YrsOut")
  aucs <- data.frame(matrix(nrow = 1, ncol = length(testYrs) + 2))
  colnames(aucs) <- c("terminalYr", "trainAUC", forecastLabels)
  aucs$terminalYr <- terminalYr
  
  # Now predict and calculate the AUC on the model training data
  train$pred <- predict(gam1, train, type = "response")
  aucs$trainAUC <- auc(train$pa, train$pred, direction = "<", quiet = TRUE)[1]
  # Predict to test data
  test$pred <- predict(gam1, test, type = "response")
  # Now loop through to get test forecast AUCs
  for (i in 1:length(testYrs)) {
    forecastYr <- testYrs[i]
    forecastData <- subset(test, year == forecastYr)
    aucs[, i + 2] <- auc(forecastData$pa, forecastData$pred, direction = "<", quiet = TRUE)[1]
  }
  # Just return AUCs for now
  out <- aucs
  return(out)
}

###########################################################################################################
# Loop through to get train and test AUCs from forecasting starting at different terminal years
# Use at least 7 years to train the model, so number of available years is
allYrs <- sort(unique(obs$year))
terminalYrs <- seq(min(allYrs) + 7, max(allYrs) - 1) # Could easily be adjusted
allAUCs <- data.frame(matrix(nrow = 0, ncol = 7))
colnames(allAUCs) <- c("terminalYr", "trainAUC", "1YrsOut", "2YrsOut", "3YrsOut", "4YrsOut", "5YrsOut")
for (j in 1:length(terminalYrs)) {
  forecastAUCs <- buildGAM(obs = obs, terminalYr = terminalYrs[j])
  allAUCs <- dplyr::bind_rows(allAUCs, forecastAUCs)
}

# Quick plot
aucsLong <- pivot_longer(allAUCs, cols = c("trainAUC", "1YrsOut", "2YrsOut", "3YrsOut", "4YrsOut", "5YrsOut"), 
                         names_to = "horizon", values_to = "auc")
# Re-order horizon as a factor
aucsLong$horizon <- factor(aucsLong$horizon, levels = c("trainAUC", "1YrsOut", "2YrsOut", "3YrsOut", "4YrsOut", "5YrsOut"))
aucsAgg <- aggregate(auc ~ horizon, aucsLong, FUN = mean, na.rm = TRUE)
aucsAgg$sd <- aggregate(auc ~ horizon, aucsLong, FUN = sd, na.rm = TRUE)[,2]
ggplot(aucsAgg) + geom_bar(aes(x = horizon, y = auc), stat = "identity") +
  geom_errorbar(aes(x = horizon, y = auc, ymin = auc - sd, ymax = auc + sd), width = 0.5) +
  scale_y_continuous(limits = c(0.5, 0.9), oob = rescale_none) + theme_bw()