###########################################################################################################
# Build an SDM (GAM or BRT) and return the model object. Just works for a binomial model at present
# k is the number of knots in the gAM, lr and tc are learning rate and tree complexity for the BRT
# Contact Barbara.Muhling@noaa.gov
###########################################################################################################

buildSDM <- function(sdmType, train, varNames, targetName, k, tc, lr) {
  if(sdmType == "gam") {
      # Build formula
      fm <- paste('s(', varNames, ', k = k', ')', sep = "", collapse = ' + ')
      fm <- as.formula(paste(targetName, '~', fm)) 
      # Train model
      mod1 <- gam(fm, data = train, family = "binomial", method = "REML", select = TRUE)
      # summary(mod1)
      # plot(mod1, scale = 0, pages = 1)
  } else if (sdmType == "brt") {
    train <- data.frame(train) # gbm.step does not like tibbles
    set.seed(1)
    mod1 <- gbm.step(data = train, gbm.x = varNames, gbm.y = targetName, 
                     tree.complexity = tc, learning.rate = lr, bag.fraction = 0.6, family = "bernoulli")
    # A simple catch for a BRT with slightly too few trees
    # This is not a substitute for providing sensible start values for tc/lr!
    # If length of training data being tested varies a lot, we will somehow need to optimize lr/tc
    ##### Flag for someone to have a go at that? #####
    if(mod1$n.trees < 1500) {
      mod1 <- gbm.step(data = train, gbm.x = varNames, gbm.y = targetName, 
                       tree.complexity = tc, learning.rate = (lr * 0.8), bag.fraction = 0.6, family = "bernoulli")
    }
  }
  
  # Return model
  return(mod1)
}