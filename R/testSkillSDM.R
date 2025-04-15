###########################################################################################################
# Test the SDM 1-year forecast skill using levels of spatial and temporal aggregation
# This function is applied to predictions to test data from a single SDM 
# Contact Barbara.Muhling@noaa.gov
###########################################################################################################

testSkillSDM <- function(mod, targetName, aucCutoff = 10) {
  # Define the test data from the list that includes the model object and test data
  test <- mod$test
  
  # Drop observations with pred = NA (from environmental data missing, e.g. outside ROMS domain)
  test <- subset(test, !is.na(pred))
  
  # Build a dataframe to take all the skill metrics: at present we're aggregating by 
  # Space: 0.1 degrees (~ 10km, native ROMS), 0.25 degrees (~ 25km), 1 degree (~ 100 km)
  # Time: monthly, seasonal
  # These could easily be adjusted 
  space = c(0.1, 0.25, 1)
  
  # Identify the months and quarters with enough observations to calculate AUC
  # and create a data.frame with all possible combinations to evaluate
  ##### (This is janky and could be improved!) #####
  months <- aggregate(eval(as.name(targetName)) ~ month + year, test, FUN = length)
  quarters <- aggregate(eval(as.name(targetName)) ~ quarter + year, test, FUN = length)
  # years <- aggregate(eval(as.name(targetName)) ~ year, test, FUN = length) # Length same as YrsToForecast
  # colnames(months)[2] <- colnames(quarters)[2] <- colnames(years)[2] <- "nObs"
  colnames(months)[3] <- colnames(quarters)[3] <- "nObs"
  months <- subset(months, nObs >= aucCutoff)
  quarters <- subset(quarters, nObs >= aucCutoff)
  months$time <- "month"
  quarters$time <- "quarter"
 
  # Expand dfs to include all options in time and space. Janky...
  months <- months[rep(1:nrow(months), length(space)),]
  months$space <- rep(space, times = (nrow(months) / length(space)))
  quarters <- quarters[rep(1:nrow(quarters), length(space)),] 
  quarters$space <- rep(space, times = (nrow(quarters) / length(space)))
  
  # Join together
  colnames(months)[1] <- colnames(quarters)[1] <- "value"
  skill <- rbind(months, quarters)
  skill$auc <- skill$nrows <- NA

  # Loop through to get AUC for each combination of spatial/temporal aggregation
  # Do this separately for multiple years contained within test (i.e. can forecast at different time horizons)
  for (i in 1:nrow(skill)) {
    spc <- skill$space[i]
    tme <- as.character(skill$time[i])
    yr <- skill$year[i]
    
    # Subset/aggregate test data accordingly. Just for one forecast year at a time
    toAgg <- subset(test, eval(as.name(skill$time[i])) == skill$value[i] & 
                          year == skill$year[i]) 
    
    # To aggregate in space, first calculate a rounding factor, then round lon/lat for aggregating
    roundFactor <- 1 / spc
    toAgg$lonrd <- round(toAgg$lon * roundFactor, 0) / roundFactor 
    toAgg$latrd <- round(toAgg$lat * roundFactor, 0) / roundFactor 
    
    # Now aggregate: let's calculate the mean of predictions, and the max of presence/absence observations
    # (that is: if any stations within spatiotemporal window have presence, then presence = 1)
    testAgg <- aggregate(eval(as.name(targetName)) ~ lonrd + latrd + eval(as.name(tme)), 
                         data = toAgg, FUN = max, na.rm = TRUE) 
    testAggPred <- aggregate(pred ~ lonrd + latrd + eval(as.name(tme)), data = toAgg, FUN = mean, na.rm = TRUE)
    # Fix colnames and join
    colnames(testAgg)[3:4] <- c(tme, targetName)
    colnames(testAggPred)[3] <- tme
    testAgg <- left_join(testAgg, testAggPred[c("lonrd", "latrd", "pred")], by = c("lonrd", "latrd"))
    
    # Calculate AUC if sufficient rows, and both presences/absence are available
    nOutcomes <- aggregate(pred ~ eval(as.name(targetName)), testAgg, FUN = length)
    if(nrow(testAgg) < aucCutoff | nrow(nOutcomes) != 2) {
      next
    } else {
      skill$auc[i] <- auc(testAgg[, 4], testAgg$pred, direction = "<", quiet = TRUE)[1]
    }
     skill$nrows[i] <- nrow(testAgg) 
  }
  # Return the df showing skill at different levels of aggregation
  return(skill)
}