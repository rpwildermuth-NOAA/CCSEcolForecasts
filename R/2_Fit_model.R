###############################################################################
#Script created by HW and modified by NLO to 1) fit and 2) Validate a BRT model
###############################################################################


#1- Select parameters
#2- Run the model 10 times so we can account for the mean and sd
#3- Validate the model (same here, better 10 iteractions so you can account for uncertainty)
#4- Historical predictions. Some people use 100% of data for predictiong, some people only the training dataset
#It dependes on how much data you have. If we have enough, I use 100%. 


################################################################################################

rm(list=ls())
#Load libraries
library(dplyr)
library(gbm)
library(dismo)
library(caret)  # For data partitioning


#read our dataset with the presences-absences and the env. data
MOM6variables=read.csv("C:\\Users\\nereo\\Dropbox (Personal)\\Nerea\\NOAA\\PROJECTS & COLLABORATIONS\\PROJECTS\\Forecast\\Swordfish_forecast\\2_build_model/swor_dataframe_MOM6_native.csv", sep=",")

data=MOM6variables
names(data)
table(data$PresAbs)
head(data)
names(data)



#filter to the roms domain before fitting the model
lat_range <- c(30, 42.5)
lon_range <- c(-130, -115)  # assuming longitudes are in degrees west (negative)

# Filter your dataframe
data <- data %>%
  filter(lat >= lat_range[1], lat <= lat_range[2],
         lon >= lon_range[1], lon <= lon_range[2])


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


#We usually run the models 10 times to have the average of the evaluation metrics and predictions in case you want to include CV
#Make function to fit 10 models, and output as a list in rds format
#Observer data only
fit.brt.n10 <- function(data, gbm.x, gbm.y, lr, iterations){
  Species_Models <- vector("list",10)
  for (i in 1:iterations){
    model <- fit.brt(data=data, gbm.x= gbm.x, gbm.y=gbm.y,lr=lr)
    Species_Models[[i]] <- model
  }
  return(Species_Models)
}



########
#Percent deviance explained: SINGLE model
#((null deviance - residual deviance)/null deviance)*100
dev_eval=function(model_object){
  null <- model_object$self.statistics$mean.null
  res <- model_object$self.statistics$mean.resid
  dev=((null - res)/null)*100 
  return(dev)
}


max(data$date)
min(data$date)
head(data)

#create a random variable and check again if the variables are falling below the random variable
#data$rand<-sample(1:100, nrow(data), replace = TRUE)

#define the env. variables names we want to use to fit the model from our dataset
gbm.x <-  c("sal", "ssh", "ssu", "ssv", "sst", "depth")  



#Fit BRT
DataInput_Fit <- data 
DataInput_Fit=na.omit(DataInput_Fit)
head(DataInput_Fit)


# Set seed for reproducibility
set.seed(123)

# Assume your dataset is named 'data' and the response variable is 'presence'
# Create an 75-25 train-test split
trainIndex <- createDataPartition(DataInput_Fit$PresAbs, p = 0.75, list = FALSE)
# Subset data
train <- DataInput_Fit[trainIndex, ]
test <- DataInput_Fit[-trainIndex, ]



Species <- "SWOR"
#quartz(width = 12, height = 8)
res1.tc3.lr01 <- fit.brt(data=train, gbm.x=gbm.x, gbm.y="PresAbs",lr=0.01) 
summary(res1.tc3.lr01) #This function give us the % importance of each env. variable to explain the distribution of our species
summary_res <- summary(res1.tc3.lr01)
setwd("C:\\Users\\nereo\\Dropbox (Personal)\\Nerea\\NOAA\\PROJECTS & COLLABORATIONS\\PROJECTS\\Forecast\\Swordfish_forecast\\3_output/native/")
write.csv(summary_res, file = "summary_results_native.csv", row.names = FALSE)

dev_eval(res1.tc3.lr01)# deviance explained by the model
summary_dev <- dev_eval(res1.tc3.lr01)
write.csv(summary_dev, file = "summary_dev_single_native.csv", row.names = FALSE)

windows(20,20)
gbm.plot(res1.tc3.lr01, smooth=TRUE, plot.layout = c(4,4), write.title=T) #response curves/ fitted functions
dev.print(png, file = "partial_curves_native.png", width = 12, height = 8, units = "in", res = 300)
dev.off()

saveRDS(res1.tc3.lr01,"C:\\Users\\nereo\\Dropbox (Personal)\\Nerea\\NOAA\\PROJECTS & COLLABORATIONS\\PROJECTS\\Forecast\\Swordfish_forecast\\3_output/native/SWOR.res1.tc3.lr03.single_native.rds")
#res1.tc3.lr01=readRDS("/Volumes/Triple_Bottom_Line/Nerea_working/Hawaii/Model_Results/Nuevo/YFT.res1.tc3.lr03.single_1res.rds")


#People usually get the normal partial curves but they are very quadratic and not smooth (like in GAMs). 
#This new function, do a boostrap and get smooth partial curves
library(ggplot2)
#setwd("/Volumes/Triple_Bottom_Line/Nerea_working/Hawaii/Model_Results/Nuevo/")
ggInfluence(res1.tc3.lr01)
dev.print(png, file = "importance_noPh.png", width = 12, height = 8, units = "in", res = 300)
dev.off()


brt1.prerun_m<- plot.gbm.4list(res1.tc3.lr01)
brt1.boot_yft <- gbm.bootstrap.functions(res1.tc3.lr01, list.predictors=brt1.prerun_m, n.reps=200)
brt1.boot_yft
saveRDS(brt1.boot_yft,"C:\\Users\\nereo\\Dropbox (Personal)\\Nerea\\NOAA\\PROJECTS & COLLABORATIONS\\PROJECTS\\Forecast\\Swordfish_forecast\\3_output/brt1.boot_yft.rds")
#


library(gridExtra)
windows(20,20)
plots=ggPD_boot_compare(gbm.object1 = res1.tc3.lr01, n.plots = 4, col.line = "darkgreen", 
                        list.4.preds1=brt1.prerun_m,  booted.preds1=brt1.boot_yft$function.preds,
                        smooth = T, col.smooth="blue",cis=c(0.025, 0.975), rug = TRUE, type.ci = "ribbon",rug.pos = "t", n.reps = 15, nrow = 2, ncol = 3)
ggsave("partial_boosted_curves.png",plot = plots, width = 40, height = 20, units = "cm")



######################################################
# Two options here: you can run the model once, and use the output for prediction
# Some people (including me( usually run 10 iteractions to account for CV and have an average performance metrics (AUCC, TSS, Sensitivity, Specificity)





#fit 10 BRT models to estimate CV
quartz(width = 12, height = 8)
res1.tc3.lr01.10models <- fit.brt.n10(data=train,gbm.x= gbm.x, gbm.y="PresAbs",lr=0.01, iterations=10)
saveRDS(res1.tc3.lr01.10models,"SWOR.res1.tc3.lr03.10iter_1res_native.rds")



#Get % deviance explained for 10 models
SDR_10models_evals <- as.data.frame(matrix(1,nrow=10,ncol=1))
counter=1
for (d in 1:10){
  eval <- (dev_eval(res1.tc3.lr01.10models[[d]]))
  SDR_10models_evals[counter,1] <- eval
  counter=counter+1
}
mean=mean(SDR_10models_evals$V1)
sd=sd(SDR_10models_evals$V1)



summary_1_mean <- mean
write.csv(summary_1_mean, file = "summary_1_mean_native.csv", row.names = FALSE)
summary_1_sd <- sd
write.csv(summary_1_sd, file = "summary_1_sd_native.csv", row.names = FALSE)


##############################################################
#3- Validation of the model
##############################################################
#We are gonna try different validations:
#Using the same dataset: Validation with full dataset (100%) and 75. vs 25 dataset (this script), and one-year out


#Leave One year Out analysis
LOO_eval <- function(DataInput, gbm.x, gbm.y, lr=lr, tc){
  DataInput <- data
  DataInput$date = as.POSIXct(DataInput$date, format = '%Y-%m-%d')
  DataInput$Year <- format(DataInput$date, "%Y")
  Evaluations_LOO <- as.data.frame(matrix(data=0,nrow=1,ncol=4))
  colnames(Evaluations_LOO) <- c("k","Deviance","AUC","TSS")
  counter=1
  for (y in min(DataInput$Year):max(DataInput$Year)){
    print(y)
    DataInput_train <- DataInput[DataInput$Year!=y,]
    DataInput_test <- DataInput[DataInput$Year==y,]
    DataInput.loo <- gbm.step(data=DataInput_train, gbm.x= gbm.x, gbm.y = c("PresAbs"), 
                              family="bernoulli", tree.complexity=tc,
                              learning.rate = lr, bag.fraction = 0.6)
    preds <- predict.gbm(DataInput.loo, DataInput_test,
                         n.trees=DataInput.loo$gbm.call$best.trees, type="response")
    dev <- calc.deviance(obs=DataInput_test$PresAbs, pred=preds, calc.mean=TRUE)
    d <- cbind(DataInput_test$PresAbs, preds)
    pres <- d[d[,1]==1,2]
    abs <- d[d[,1]==0,2]
    if(length(pres)>0 & length(abs)>0){
      e <- evaluate(p=pres, a=abs)
      
      Evaluations_LOO[counter,1] <- y
      Evaluations_LOO[counter,2] <- dev
      Evaluations_LOO[counter,3] <- e@auc
      Evaluations_LOO[counter,4] <- max(e@TPR + e@TNR-1)
      counter=counter+1 
    }
  }
  return(Evaluations_LOO)}


#######
#Make function to 75/25 split AUC test. 
#This is to 'repeat' what Elliot thinks they did for EcoCast (Kylie disagrees, see next function)
eval_7525 <- function(DataInput, gbm.x, gbm.y, lr, tc=tc, family){
  DataInput <- DataInput_Fit
  Evaluations_7525 <- as.data.frame(matrix(data=0,nrow=1,ncol=5))
  colnames(Evaluations_7525) <- c("Deviance","AUC","TSS","Sensitivity", "Specificity")
  DataInput_bound <- floor((nrow(DataInput)/4)*3)         #define % of training and test set
  DataInput_train<- DataInput[sample(nrow(DataInput),DataInput_bound),]
  DataInput_test<- sqldf('SELECT * FROM DataInput EXCEPT SELECT * FROM DataInput_train')
  dim(DataInput_test)
  dim(DataInput_train)
  DataInput.kfolds <- gbm.step(data=DataInput_train, gbm.x= gbm.x, gbm.y = "PresAbs", 
                               family="bernoulli", tree.complexity=tc,
                               learning.rate = lr, bag.fraction = 0.6)
  preds <- predict.gbm(DataInput.kfolds, DataInput_test,
                       n.trees=DataInput.kfolds$gbm.call$best.trees, type="response")
  dev <- calc.deviance(obs=DataInput_test$PresAbs, pred=preds, calc.mean=TRUE)
  d <- cbind(DataInput_test$PresAbs, preds)
  pres <- d[d[,1]==1,2]
  abs <- d[d[,1]==0,2]
  e <- evaluate(p=pres, a=abs)
  Evaluations_7525[1,1] <- dev
  Evaluations_7525[1,2] <- e@auc
  Evaluations_7525[1,3] <- max(e@TPR + e@TNR-1)
  Evaluations_7525[1,4] <- mean(e@TPR)
  Evaluations_7525[1,5] <- mean(e@TNR)
  return(Evaluations_7525)}


#Make function to 100/100  AUC test
eval_100_percent <- function(dataInput, gbm.x, gbm.y, lr=lr, tc){
  DataInput <- dataInput
  Evaluations_100_percent <- as.data.frame(matrix(data=0,nrow=1,ncol=3))
  colnames(Evaluations_100_percent) <- c("Deviance","AUC","TSS")
  DataInput_train<- DataInput
  DataInput_test<- DataInput
  DataInput.kfolds <- gbm.step(data=DataInput_train, gbm.x= gbm.x, gbm.y = gbm.y, 
                               family="bernoulli", tree.complexity=tc,
                               learning.rate = lr, bag.fraction = 0.6)
  preds <- predict.gbm(DataInput.kfolds, DataInput_test,
                       n.trees=DataInput.kfolds$gbm.call$best.trees, type="response")
  dev <- calc.deviance(obs=DataInput_test$PresAbs, pred=preds, calc.mean=TRUE)
  d <- cbind(DataInput_test$PresAbs, preds)
  pres <- d[d[,1]==1,2]
  abs <- d[d[,1]==0,2]
  e <- evaluate(p=pres, a=abs)
  Evaluations_100_percent[1,1] <- dev
  Evaluations_100_percent[1,2] <- e@auc
  Evaluations_100_percent[1,3] <- max(e@TPR + e@TNR-1)
  
  return(Evaluations_100_percent)}

library(sqldf)
SDR.loo.eval <-LOO_eval(DataInput_Fit,gbm.x=gbm.x, gbm.y="PresAbs",lr=0.01, tc=3)
SDR.7525.eval <- eval_7525(DataInput_Fit,gbm.x=gbm.x, gbm.y="PresAbs",lr=0.01, tc=3)
SDR.100.eval <- eval_100_percent(DataInput_Fit,gbm.x=gbm.x, gbm.y="PresAbs",lr=0.01, tc=3)



saveRDS(SDR.loo.eval,paste("SWOR_SDR.loo.eval_native.rds"))
saveRDS(SDR.7525.eval,paste("SWOR_SDR.7525.eval_native.rds"))
saveRDS(SDR.100.eval,paste("SWOR_SDR.100.eval_native.rds"))


SDR.7525.eval <- SDR.7525.eval
write.csv(SDR.7525.eval, file = "SDR.7525.eval_native.csv", row.names = FALSE)

SDR.loo.eval <- SDR.loo.eval
write.csv(SDR.loo.eval, file = "SDR.loo.eval_native.csv", row.names = FALSE)

SDR.100.eval <- SDR.100.eval
write.csv(SDR.100.eval, file = "SDR.100.eval_native.csv", row.names = FALSE)





