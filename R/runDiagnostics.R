###########################################################################################################
# A quick function to return some diagnostics from a fitted SDM
# More metrics could be added as desired
# Currently works for a GAM or BRT (binomial only)
# Contact bmuhling@ucsc.edu
###########################################################################################################

runDiagnostics <- function(mod, max.trees) {
  modelObject <- mod$sdm
  testData <- mod$test
  
  if("gam" %in% class(modelObject)) {
    # print("model is a GAM")
    trainingData <- modelObject$model
    # Did model converge?
    cnv <- modelObject$converged # TRUE or FALSE
    # Calculate the training data AUC, no CV for GAMs at the moment
    trainingData$pred <- predict(modelObject, trainingData, type = "response")
    trainAUC <- auc(trainingData[, 1], trainingData$pred, direction = "<", quiet = TRUE)[1] 
    cvAUC <- NA
  }
  
  if("gbm" %in% class(modelObject)) {
    # print("model is a BRT")
    trainingData <- modelObject$gbm.call$dataframe 
    # Did model converge? 
    # BRTs pretty much always converge, but am testing here whether n.trees = max.trees, which indicates that
    # lr is too low and/or tc is too high. Will throw warning when running gbm.step:
    # "maximum tree limit reached - results may not be optimal: refit with faster learning rate or increase maximum number of trees" 
    # but warning not saved to model object, annoyingly neither is max.trees, so need to pass it to this fn
    nTrees <- modelObject$n.trees # Final fitted trees 
    cnv <- ifelse(nTrees == max.trees, FALSE, TRUE)
    # Extract the training data AUC and cross-validation AUC for BRT
    trainAUC <- modelObject$self.statistics$discrimination
    cvAUC <- modelObject$cv.statistics$discrimination.mean
  }
  # Combine and return
  evalOutputs <- list("convergence" = cnv, "trainingAUC" = trainAUC, "crossValidationAUC" = cvAUC)
  return(evalOutputs)
}