######SCRATCH TO READ THE SWORDFISH MODEL AND TRY TO APPLY BARB'S FUNCTION FOR FORECAST TESTING 
###########################################################################################################
# Evaluate the fit of an SDM on withheld test data
# Testing withholding the last 1 - 10 years of available time-series
# This example shows a BRT for swordfish presence/absence
# Contact nerea.lezama-ochoa@noaa
###########################################################################################################

#Load libraries
rm(list = ls())
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

brt_model <- readRDS("C:\\Users\\nereo\\Dropbox (Personal)\\Nerea\\NOAA\\PROJECTS & COLLABORATIONS\\PROJECTS\\Forecast\\Swordfish_forecast\\3_output/native/SWOR.res1.tc3.lr03.single_native.rds")

dataframe=brt_model[["gbm.call"]][["dataframe"]]
#write.table(dataframe, file="C:\\Users\\nereo\\OneDrive\\Escritorio\\Barb/swor_dataframe.csv", sep=",")
head(dataframe)
max(dataframe$date)
min(dataframe$date)
# If observations don't have a "year" column, add one
dataframe$year <- year(dataframe$date)
# I'm missing predictor variables for 2024, so trimming that year off
dataframe <- subset(dataframe, year <= 2017)
# Define the target variable
dataframe$pa <- dataframe$PresAbs


#Function to fit brts
fit.brt <- function(data, gbm.x, gbm.y, lr=lr){
  y <- gbm.step(data=data, 
                gbm.x = gbm.x,     ### indices of columns holding predictor variables, as vector (e.g. c(8,10:13))
                gbm.y = "PresAbs",     ### index of column holding response variable (e.g. presabs)
                family = "bernoulli",
                tree.complexity = 3, ### this is the complexity of the interactions that the model will fit.  We should restrict <=3
                learning.rate = lr,  ### this needs to be optimised to end up with >1000 trees, from the approximate range 0.01-0.05.  
                bag.fraction = 0.60)  ### recommended by Elith, amount of input data used each time
  y$gbm.call$dataframe<-data
  return(y)}

gbm.x <-  c("sal", "ssh", "ssu", "ssv", "sst","depth")  


head(dataframe)

buildBRT <- function(dataframe, terminalYr, k) {
  # Define the training and test forecast year
  train <- subset(dataframe, year <= terminalYr)
  test <- subset(dataframe, year > terminalYr & year <= (terminalYr + 5)) # May be less than 5 years if missing years, e.g. 2020
  
  # Now build a BRT 
  brt1 <- fit.brt(data=train, gbm.x=gbm.x, gbm.y="pa",lr=0.01)
  # summary(brt1)
  
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
  train$pred <- predict(brt1, train, type = "response")
  aucs$trainAUC <- auc(train$pa, train$pred, direction = "<", quiet = TRUE)[1]
  # Predict to test data
  test$pred <- predict(brt1, test, type = "response")
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
allYrs <- sort(unique(dataframe$year))
terminalYrs <- seq(min(allYrs) + 7, max(allYrs) - 1) # Could easily be adjusted
allAUCs <- data.frame(matrix(nrow = 0, ncol = 7))
colnames(allAUCs) <- c("terminalYr", "trainAUC", "1YrsOut", "2YrsOut", "3YrsOut", "4YrsOut", "5YrsOut")
for (j in 1:length(terminalYrs)) {
  forecastAUCs <- buildBRT(dataframe = dataframe, terminalYr = terminalYrs[j])
  allAUCs <- dplyr::bind_rows(allAUCs, forecastAUCs)
}






# Quick plot
aucsLong <- pivot_longer(allAUCs, cols = c("trainAUC", "1YrsOut", "2YrsOut", "3YrsOut", "4YrsOut", "5YrsOut"), 
                         names_to = "horizon", values_to = "auc")
# Re-order horizon as a factor
aucsLong$horizon <- factor(aucsLong$horizon, levels = c("trainAUC", "1YrsOut", "2YrsOut", "3YrsOut", "4YrsOut", "5YrsOut"))
aucsAgg <- aggregate(auc ~ horizon, aucsLong, FUN = mean, na.rm = TRUE)
aucsAgg$sd <- aggregate(auc ~ horizon, aucsLong, FUN = sd, na.rm = TRUE)[,2]
windows(20,20)
ggplot(aucsAgg) + geom_bar(aes(x = horizon, y = auc), stat = "identity") +
  geom_errorbar(aes(x = horizon, y = auc, ymin = auc - sd, ymax = auc + sd), width = 0.5) +
  scale_y_continuous(limits = c(0.5, 0.9), oob = rescale_none) + theme_bw()


setwd("C:\\Users\\nereo\\Dropbox (Personal)\\Nerea\\NOAA\\PROJECTS & COLLABORATIONS\\PROJECTS\\Forecast\\Swordfish_forecast\\5_testing_skil\\")
write.csv(aucsAgg, file = "auc_native.csv", row.names = FALSE)




