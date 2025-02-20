###########################################################################################################
# Build an SDM (GAM or BRT) by calling buildSDM, return the model object and predictions to withheld
# test data (1 year forecast)
# Contact Barbara.Muhling@noaa.gov
###########################################################################################################

scoreSDM <- function(subObs, sdmType, varNames, targetName, k, tc, lr) {
  # Define the training and test forecast years
  yrs <- unique(sort(subObs$year)) # Years in the observational dataset
  terminalYr <- max(yrs) - 1 # The last year of training data
  train <- subset(obs, year <= terminalYr)
  test <- subset(obs, year == (terminalYr + 1)) # Observations the year after the training data ends
  
  # Sometimes a whole year of observations is missing (e.g. 2020), or there are very few observations
  # In that case, stop and return NA
  if(nrow(test) < 10) { # Could use another cutoff, here it's 10 observations
    return(NA)
  }
  
  # Otherwise, build an SDM using helper function
  source("./R/buildSDM.R")
  mod1 <- buildSDM(sdmType = sdmType, train = train, varNames = varNames, targetName = targetName, k = k, tc = tc, lr = lr)
  # summary(mod1) # If you want to check convergence etc. But GAMs/BRTs nearly always converge unless parameters v inappropriate
  # gbm.step prints model convergence progress as it goes, so you'll see if the number of trees is too small (< ~ 1500)
  
  # Score the test dataset: the one year of data following the training data (1 year forecast)
  test$pred <- predict(mod1, test, type = "response")
  
  # For now, returning model object (which contains the training data), as well as the test/forecast data
  # For a large training dataset, this could result in a very large object though
  out <- list("sdm" = mod1, "test" = test) 
  return(out)
}