

# # Script to extract the MOM6 variables in native grid from my swordfish dataset
#Script created by Nerea Lezama-Ochoa (nerea.lezama-ochoa@noaa.gov)

rm(list=ls())
require(lubridate)
require(ncdf4)
require(dplyr)
library(sp)


outputs="C:\\Users\\nereo\\Dropbox (Personal)\\Nerea\\NOAA\\PROJECTS & COLLABORATIONS\\PROJECTS\\Forecast\\Swordfish_forecast\\3_output/"

#read dataframe for swordfish
swor=read.csv("C:\\Users\\nereo\\Dropbox (Personal)\\Nerea\\NOAA\\PROJECTS & COLLABORATIONS\\PROJECTS\\Forecast\\Swordfish_forecast\\1_extract_variables/swor_dataframe.csv", sep=",")
head(swor)
swor <- swor[, c("lon", "lat", setdiff(names(swor), c("lon", "lat")))]
swor <- swor[, c("lon", "lat", "dt", "PresAbs")]
swor$dt <- as.Date(swor$dt, format = "%m/%d/%Y")
swor$date=swor$dt
max(swor$date)
min(swor$date)
head(swor)

#filter to the dates we have
swor <- swor %>%
  filter(dt > as.Date("1993-01-01") & dt < as.Date("2017-01-31"))


swor <- swor[, -3]
dim(swor)
head(swor)



desired.diameter <- 0.1 
timestep <- "daily"
func = "mean" # "mean" or "sd"
nc.path <- "C:\\Users\\nereo\\Dropbox (Personal)\\Nerea\\NOAA\\PROJECTS & COLLABORATIONS\\PROJECTS\\Forecast\\Swordfish_forecast\\1_extract_variables\\MOM6/"



#salinity
varName <- "sos"
desired.diameter <- 0.1 
timestep <- "daily"
func = "mean" # "mean" or "sd"

sal <- getMOM6_raw(points = swor, varName = varName, desired.diameter = desired.diameter, func = func,
               timestep = timestep, nc.path = nc.path)

#sal <- sal %>% select(-lon, -lat, -date)
head(sal)


#ssh
varName <- "ssh"
desired.diameter <- 0.1 
timestep <- "daily"
func = "mean" # "mean" or "sd"

ssh <- getMOM6_raw(points = swor, varName = varName, desired.diameter = desired.diameter, func = func,
               timestep = timestep, nc.path = nc.path)

ssh <- ssh[, !(names(ssh) %in% c("lon", "lat", "date", "PresAbs"))]
head(ssh)


#ssu
varName <- "ssu"
ssu <- getMOM6_raw(points = swor, varName = varName, desired.diameter = desired.diameter, func = func,
               timestep = timestep, nc.path = nc.path)


ssu <- ssu[, !(names(ssu) %in% c("lon", "lat", "date", "PresAbs"))]
head(ssu)


#ssv
varName <- "ssv"
desired.diameter <- 0.1 
timestep <- "daily"
func = "mean" # "mean" or "sd"
nc.path <- "C:\\Users\\nereo\\Dropbox (Personal)\\Nerea\\NOAA\\PROJECTS & COLLABORATIONS\\PROJECTS\\Forecast\\Swordfish_forecast\\1_extract_variables\\MOM6//"

ssv <- getMOM6_raw(points = swor, varName = varName, desired.diameter = desired.diameter, func = func,
                      timestep = timestep, nc.path = nc.path)

ssv <- ssv[, !(names(ssv) %in% c("lon", "lat", "date", "PresAbs"))]
head(ssv)


#sst
varName <- "tos"
desired.diameter <- 0.1 
timestep <- "daily"
func = "mean" # "mean" or "sd"

sst <- getMOM6_raw(points = swor, varName = varName, desired.diameter = desired.diameter, func = func,
                      timestep = timestep, nc.path = nc.path)

sst <- sst[, !(names(sst) %in% c("lon", "lat", "date", "PresAbs"))]
head(sst)




#O2
#varName <- "btm_o2"
#desired.diameter <- 0.1 # Note! This means a 0.7x0.7 degree box
#timestep <- "daily"
#func = "mean" # "mean" or "sd"

#O2 <- getMOM6(points = swor, varName = varName, desired.diameter = desired.diameter, func = func,
#             timestep = timestep, nc.path = nc.path)

#O2 <- O2[, !(names(O2) %in% c("lon", "lat", "date", "PresAbs"))]
#head(O2)


#depth
varName <- "deptho"
desired.diameter <- 0.1 
#timestep <- "daily"
func = "mean" # "mean" or "sd"
nc.file = "C:\\Users\\nereo\\Dropbox (Personal)\\Nerea\\NOAA\\PROJECTS & COLLABORATIONS\\PROJECTS\\Forecast\\Swordfish_forecast\\1_extract_variables\\MOM6/daily/ocean_daily.19930101-20191231.deptho.nc"

depth <- getMOM6_static_depth(points = swor, varName = varName, nc.file = nc.file, desired.diameter = desired.diameter, fun=fun)
head(depth)





#merge all the variables
final=cbind(sal, ssh, ssu, ssv, sst, depth)
head(final)


#select the columns
MOM6variables <- final[, c("lon", "lat", "PresAbs", "date", 
                           "sos_mean_0.1_daily", "ssh_mean_0.1_daily", 
                           "ssu_mean_0.1_daily", "ssv_mean_0.1_daily", 
                           "tos_mean_0.1_daily", "deptho_value")]
head(MOM6variables)

#change the column names to make it easier
names(MOM6variables) <- c("lon", "lat", "PresAbs", "date",  "sal", "ssh", "ssu", "ssv", "sst", "depth")  

#save the dataset
write.table(MOM6variables, file="C:\\Users\\nereo\\Dropbox (Personal)\\Nerea\\NOAA\\PROJECTS & COLLABORATIONS\\PROJECTS\\Forecast\\Swordfish_forecast\\2_build_model/swor_dataframe_MOM6_native.csv", sep=",")


 