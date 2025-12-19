# Code from Nerea Ochoa Lazama to evaluate SDM sensitivity to data inputs using
# either leave-one-out cross-validation, or k-fold cross-validation

#!!RW: Currently only defined for BRT models


#Leave One year Out analysis
LOO_eval <- function(DataInput, gbm.x, gbm.y, lr, tc){
  DataInput$date <- as.POSIXct(DataInput$date, format = '%Y-%m-%d')
  DataInput$Year <- format(DataInput$date, "%Y")
  Evaluations_LOO <- as.data.frame(matrix(data=0,nrow=1,ncol=4))
  colnames(Evaluations_LOO) <- c("k","Deviance","AUC","TSS")
  counter=1
  for (y in min(DataInput$Year):max(DataInput$Year)){
    print(y)
    DataInput_train <- DataInput[DataInput$Year!=y,]
    DataInput_test <- DataInput[DataInput$Year==y,]
    DataInput.loo <- gbm.step(data=DataInput_train, gbm.x= gbm.x, gbm.y = gbm.y, 
                              family="bernoulli", tree.complexity=tc,
                              learning.rate = lr, bag.fraction = 0.6)
    preds <- predict.gbm(DataInput.loo, DataInput_test,
                         n.trees=DataInput.loo$gbm.call$best.trees, type="response")
    # !!RW: Is this the same function as the dev_eval() defined in FitModel?
    dev <- calc.deviance(obs=DataInput_test[, gbm.y], pred=preds, calc.mean=TRUE)
    d <- cbind(DataInput_test[, gbm.y], preds)
    pres <- d[d[,1]==1,2]
    abs <- d[d[,1]==0,2]
    if(length(pres)>0 & length(abs)>0){
      e <- dismo::evaluate(p=pres, a=abs) #!!RW: dismo package?
      
      Evaluations_LOO[counter,1] <- y
      Evaluations_LOO[counter,2] <- dev
      Evaluations_LOO[counter,3] <- e@auc
      Evaluations_LOO[counter,4] <- max(e@TPR + e@TNR-1)
      counter=counter+1 
    }
  }
  return(Evaluations_LOO)
  }


#######
#Make function to 75/25 split AUC test. 
#This is to 'repeat' what Elliot thinks they did for EcoCast (Kylie disagrees, see next function)
eval_7525 <- function(DataInput, gbm.x, gbm.y, lr, tc=tc, family, retainPct = 0.75){
  Evaluations_7525 <- as.data.frame(matrix(data=0,nrow=1,ncol=5))
  colnames(Evaluations_7525) <- c("Deviance","AUC","TSS","Sensitivity", "Specificity")
  
  trainIndex <- createDataPartition(DataInput[, gbm.y], p = retainPct, list = FALSE)
  # Subset data
  DataInput_train <- DataInput[trainIndex, ]
  print(dim(DataInput_train))
  if(retainPct == 1){
    DataInput_test <- DataInput_train
  } else {
    DataInput_test <- DataInput[-trainIndex, ]
  }
 
  DataInput.kfolds <- gbm.step(data=DataInput_train, gbm.x= gbm.x, gbm.y = gbm.y, 
                               family="bernoulli", tree.complexity=tc,
                               learning.rate = lr, bag.fraction = 0.6)
  preds <- predict.gbm(DataInput.kfolds, DataInput_test,
                       n.trees=DataInput.kfolds$gbm.call$best.trees, type="response")
  dev <- calc.deviance(obs=DataInput_test[, gbm.y], pred=preds, calc.mean=TRUE)
  d <- cbind(DataInput_test[, gbm.y], preds)
  pres <- d[d[,1]==1,2]
  abs <- d[d[,1]==0,2]
  e <- dismo::evaluate(p=pres, a=abs)
  Evaluations_7525[1,1] <- dev
  Evaluations_7525[1,2] <- e@auc
  Evaluations_7525[1,3] <- max(e@TPR + e@TNR-1)
  Evaluations_7525[1,4] <- mean(e@TPR)
  Evaluations_7525[1,5] <- mean(e@TNR)
  return(Evaluations_7525)
  }
