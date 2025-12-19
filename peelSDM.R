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
  #evalOutputs <- runDiagnostics(mod = mod, max.trees = max.trees)
  
  # Return results: at present, returning as a list
  if(includePersistence == TRUE) {
    out <- list("sdm" = mod$sdm, "test" = mod$test, #"modelEval" = evalOutputs, 
                "sdmSkill" = sdmSkill, "sdmSkillPers" = sdmSkillPers)
  } else {
    out <- list("sdm" = mod$sdm, "test" = mod$test, #"modelEval" = evalOutputs, 
                "sdmSkill" = sdmSkill)
  }
  return(out)
}
