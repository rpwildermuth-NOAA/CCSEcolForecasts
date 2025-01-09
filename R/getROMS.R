###########################################################################################################
# Overhauled script for extracting ROMS environmental variable to provided points
# Basically combines older "getROMS and "XtractROMS" scritps together, and runs for a selected variable
# and function (mean, sd), instead of always running everything together
# Based on previous code from Karin Forney, Heather Welch, Steph Brodie, and others
# Requires an input data frame with columns labeled "lon" (degrees west), "lat", and "date" (formatted as dates)
# Note! In response to past inconsistencies, desired.resolution has been replaced with desired.diameter
# Thus desired.diameter = 0.7 (degrees) results in a 7x7 (= 49 pixel) box around the central location
# Contact Barbara.Muhling@noaa.gov
###########################################################################################################

require(ncdf4)
require(lubridate)

###########################################################################################################
# # Testing function
# points <- data.frame(expand.grid("lon" = seq(-132, -122, by = 2), "lat" = seq(32, 40, by = 2),
#                                  "date" = seq(as.Date("2009-06-15"), as.Date("2017-06-15"), by = "month")))
# varName <- "sst"
# desired.diameter <- 0.7 # Note! This means a 7x7 or 49 pixel box
# func <- "sd" # mean or sd
# histPath <- "F:/roms/hist" # As some folks might have them mixed together, I use subfolders...
# nrtPath <- "F:/roms/nrtComplete"
# testroms <- getROMS(points = points, varName = varName, desired.diameter = desired.diameter, func = func, 
#                 histPath = histPath, nrtPath = nrtPath)

###########################################################################################################
getROMS <- function(points, varName, desired.diameter, func = "mean", histPath, nrtPath) {
  # Check that the points have columns called "lon", "lat", and "date"
  if(!"lon" %in% colnames(points) | !"lat" %in% colnames(points) | !"date" %in% colnames(points)) {
    stop("Input points data frame needs columns named lon, lat, and date")
  }
  # Check that date column is formatted as such
  if(!inherits(points$date, 'Date')) {
   stop("Date column is not formatted as dates")
  }
  
  #########################################################################################################
  ################################ Identify the correct netcdfs ###########################################
  #########################################################################################################
  # List the correct ROMS netcdfs: one each for the historical (1980-2010) and near-real-time (2011-present) time periods
  # Should work whether stored in same or different directories, 
  # as long as e.g. several versions of NRT netcdfs aren't mixed in together
  histFile <- intersect(list.files(path = histPath, pattern = "1980.*nc"),
                        list.files(path = histPath, pattern = varName))
  nrtFile <- intersect(list.files(path = nrtPath, pattern = "2011.*nc"),
                       list.files(path = nrtPath, pattern = varName))
  
  # Catch for BV/BF being called different things: in my hist ROMs filename contains "BV", NRT contains "bbv_200"
  if(varName == "bv" | varName == "BV" | varName == "bv_200" | varName == "bf" | varName == "BV_frequency") {
    histFile <- intersect(list.files(path = histPath, pattern = "1980.*nc"),
                          list.files(path = histPath, pattern = "bv|BV|bv_200|bf"))
    nrtFile <- intersect(list.files(path = nrtPath, pattern = "2011.*nc"),
                          list.files(path = nrtPath,  pattern = "bv|BV|bv_200|bf"))
    # And specify the correct name for the variable!
    varNameHist <- "BV_frequency"
    varNameNrt <- "bbv_200"
  } else {
    varNameHist <-varNameNrt <- varName
  }
  
  # Print warnings/stop and show error if there's not 1 file each
  if(length(c(histFile, nrtFile)) != 2) {
    stop("Incorrect number of files: there should be one historical and one NRT file")
  }

  #########################################################################################################
  #################################### Open the correct netcdfs ###########################################
  #########################################################################################################
  # Open netcdfs required based on date range of input points
  minDatePnts <- min(points$date)
  maxDatePnts <- max(points$date)
  if(minDatePnts <= as.Date("2010-12-31")) {
    histnc <- nc_open(paste0(histPath, "/", histFile))
    histOpen <- TRUE
  }
  if(maxDatePnts >= as.Date("2011-01-01")) {
    nrtnc <- nc_open(paste0(nrtPath, "/", nrtFile))
    nrtOpen <- TRUE
  }
  
  # Define the lon, lat, and time dimensions. This is fiddly as some files have different dimnames
  getDims <- function(nc) {
    nm <- names(nc$var)
    if("lon_rho" %in% nm) {
      lat <- ncvar_get(nc, "lat_rho"); lat <- lat[1,]
      lon <- ncvar_get(nc, "lon_rho"); lon <- lon[,1]  
    } else {
      lat <- ncvar_get(nc, "lat"); lat <- lat[1,]
      lon <- ncvar_get(nc, "lon"); lon <- lon[,1]
    }
    # Some newer files just have time, not date
    if("year" %in% nm) {
      yr <- ncvar_get(nc, "year")
      mth <- ncvar_get(nc, "month")
      day <- ncvar_get(nc, "day")
      time <- as.Date(as.POSIXct(paste(yr, mth, day, sep = "-"), tz = "UTC"))
    } else {
      tim <- ncvar_get(nc, "time")
      tim2 <- as.POSIXct((tim / 86400) + as.Date("2011-01-02"))
      yr <- year(tim2)
      mth <- month(tim2)
      day <- day(tim2)
      time <- as.Date(as.POSIXct(paste(yr, mth, day, sep = "-"), tz = "UTC"))
    }
    return(out <- list("lon" = lon, "lat" = lat, "time" = time))
  }
  
  # Run fn for whichever of hist and nrt netcdfs you opened
  if(histOpen == TRUE) {
    histDims <- getDims(nc = histnc)
  }
  if(nrtOpen == TRUE) {
    nrtDims <- getDims(nc = nrtnc)
  }
  
  #########################################################################################################
  ################################# Extract from netcdfs at points ########################################
  #########################################################################################################
  # Loop through points to extract desired variable at correct spatial resolution
  # Also get mean or stdev, depending on inputs to fn
  # Add new columns as needed to receive calculated values
  # The number in the column title is the diameter in degrees
  if(func == "mean") {
    points$mn <- NA
    colnames(points)[ncol(points)] <- paste0(varName, "_", "mean", "_", desired.diameter)
    points$count <- NA
    colnames(points)[ncol(points)] <- paste0(varName, "_", "count", "_", desired.diameter)
  }
  if(func == "sd") {
    points$sd <- NA
    colnames(points)[ncol(points)] <- paste0(varName, "_", "sd", "_", desired.diameter)
    points$count <- NA
    colnames(points)[ncol(points)] <- paste0(varName, "_", "count", "_", desired.diameter)
  }
  
  # Check that point is within spatiotemporal range (including pixel window)
  minLon <- -134 + ((desired.diameter - 0.1) / 2) 
  maxLon <- -115 - ((desired.diameter - 0.1) / 2)
  minLat <- 30 + ((desired.diameter - 0.1) / 2)
  maxLat <- 48 - ((desired.diameter - 0.1) / 2)
  minDate <- as.Date("1980-01-11")
  maxDate <- dplyr::if_else(nrtOpen == TRUE, max(nrtDims$time), as.Date("2010-12-31"))
  # Loop through points
  for (i in 1:nrow(points)) {
    if(points$lon[i] < minLon | points$lon[i] > maxLon | points$lat[i] < minLat | points$lat[i] > maxLat |
        points$date[i] < minDate | points$date[i] > maxDate) {
      next
    }
    # If point is in spatiotemporal range, proceed
    # If point is in historical date range:
    if(points$date[i] %in% histDims$time) {
      xdate <- which(points$date[i] == histDims$time)
      c <- which.min(abs(histDims$lon - points$lon[i]))
      c_low <- which.min(abs(histDims$lon - (points$lon[i] - (desired.diameter - 0.1) / 2))) 
      c_up <- which.min(abs(histDims$lon - (points$lon[i] + (desired.diameter - 0.1) / 2))) 
      r <- which.min(abs(histDims$lat - points$lat[i]))
      r_low <- which.min(abs(histDims$lat - (points$lat[i] - (desired.diameter - 0.1) / 2)))
      r_up <- which.min(abs(histDims$lat - (points$lat[i] + (desired.diameter - 0.1) / 2)))
      numcols = abs(c_up - c_low) + 1 
      numrows = abs(r_up - r_low) + 1
      
      # Extract
      data.var  <-  ncvar_get(histnc, varNameHist, start = c(c_low, r_low, xdate), # Point or 2D matrix
                                count = c(numcols, numrows, 1), verbose = FALSE)
      if(!is.na(mean(data.var, na.rm = TRUE))) {
        if (func == "mean") {
          points[i, ncol(points) - 1] <- mean(data.var, na.rm = TRUE)
          points[i, ncol(points)] <- sum(!is.na(data.var)) 
        } else if (func == "sd") {
          points[i, ncol(points) - 1] <- sd(data.var, na.rm = TRUE)
          points[i, ncol(points)] <- sum(!is.na(data.var))
        }
      }
    } 
      
    # If point is in NRT date range: 
     if(points$date[i] %in% nrtDims$time) {
       xdate <- which(points$date[i] == nrtDims$time)
       c <- which.min(abs(nrtDims$lon - points$lon[i]))
       c_low <- which.min(abs(nrtDims$lon - (points$lon[i] - (desired.diameter - 0.1) / 2))) 
       c_up <- which.min(abs(nrtDims$lon - (points$lon[i] + (desired.diameter - 0.1) / 2))) 
       r <- which.min(abs(nrtDims$lat - points$lat[i]))
       r_low <- which.min(abs(nrtDims$lat - (points$lat[i] - (desired.diameter - 0.1) / 2)))
       r_up <- which.min(abs(nrtDims$lat - (points$lat[i] + (desired.diameter - 0.1) / 2)))
       numcols = abs(c_up - c_low) + 1 
       numrows = abs(r_up - r_low) + 1
       
       # Extract
       data.var  <-  ncvar_get(nrtnc, varNameNrt, start = c(c_low, r_low, xdate), # Point or 2D matrix
                               count = c(numcols, numrows, 1), verbose = FALSE)
       if(!is.na(mean(data.var, na.rm = TRUE))) {
         if (func == "mean") {
           points[i, ncol(points) - 1] <- mean(data.var, na.rm = TRUE)
           points[i, ncol(points)] <- sum(!is.na(data.var)) 
         } else if (func == "sd") {
           points[i, ncol(points) - 1] <- sd(data.var, na.rm = TRUE)
           points[i, ncol(points)] <- sum(!is.na(data.var)) 
         }
       }
     } 
    # Status indicator for large datasets
    if((i - 100) %% 100 == 0) {
      print(paste0(i, " points complete"))
    }
  } 
  # Close open files
  if(histOpen == TRUE) {
      nc_close(histnc)
  }
  if(nrtOpen == TRUE) {
      nc_close(nrtnc) 
  }
  return(points)
} 
