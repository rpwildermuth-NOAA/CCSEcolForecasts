#####################################################################################
# Extract CMEMS level 4 daily chl interpolated product
# Requires global 4km netcdfs downloaded from 
# https://data.marine.copernicus.eu/product/OCEANCOLOUR_GLO_BGC_L4_MY_009_104/description
# and stored locally (nc.path)
# Requires an input data frame with columns labeled "lon" (degrees west), "lat", and "date" (formatted as dates)
# Note! In response to past inconsistencies, desired.resolution has been replaced with desired.diameter
# Thus desired.diameter = 0.7 results in a 0.7 x 0.7 degree box around the central location
# Unlike ROMS, CMEMS chl resolution is defined in km, not degrees. So pixel sizes will vary slightly with latitude 
# To keep workflow consistent with existing FRD/ESD workflows, we assume that pixel sizes are 0.0416667 X 0.0416667
# This may not be accurate at high latitudes!
# Contact Barbara.Muhling@noaa.gov
#####################################################################################

require(lubridate)
require(ncdf4)

###########################################################################################################
# # Testing function
# points <- data.frame("lon" = c(-130, -125, -124, -130, -125, -124, -130, -125, -124, -180),
#                      "lat" = c(45, 34, 32, 40, 45, 34, 32, 40, 36, 41),
#                      "date" = seq(as.Date("2008-08-15"), as.Date("2017-08-15"), by = "year"))
# desired.diameter <- 0.7 # Note! This means a 0.7x0.7 degree box
# func <- "mean" # mean or sd
# nc.path <- "F:/"

###########################################################################################################

getCMEMS_l4chl <- function(points, desired.diameter, func = "mean", nc.path) {
  # Check that the points have columns called "lon", "lat", and "date"
  if(!"lon" %in% colnames(points) | !"lat" %in% colnames(points) | !"date" %in% colnames(points)) {
    stop("Input points data frame needs columns named lon, lat, and date")
  }
  # Check that date column is formatted as such
  if(!inherits(points$date, 'Date')) {
    stop("Date column is not formatted as dates")
  }
  
  # If desired.diameter is smaller than native resolution, set to native resolution
  if(desired.diameter < 0.0417) {
    desired.diameter <- 0.0417
  }
  
  # Define point locations, add a lon360 column (degrees east)
  # This allows us to extract points across the dateline if needed
  points$lon360 <- ifelse(points$lon < 0, points$lon + 360, points$lon)
  
  # Define date variables (required to open correct netcdf)
  fishdate <- as.Date(points$date)
  fishyear <- year(fishdate)
  fishday  <- formatC(day(fishdate), width = 2, flag = 0)
  fishmonth <- formatC(month(fishdate), width = 2, flag = 0)

  # Loop through points to extract desired variable at correct spatial resolution
  # Also get mean or stdev, depending on inputs to fn
  # Add new columns as needed to receive calculated values
  # The number in the column title is the diameter in degrees
  if(func == "mean") {
    points$mn <- NA
    colnames(points)[ncol(points)] <- paste0("chl_mean", "_", desired.diameter)
  }
  if(func == "sd") {
    points$sd <- NA
    colnames(points)[ncol(points)] <- paste0("chl_sd", "_", desired.diameter)
  }
  
  # Now extract for all data points
  for(i in 1:nrow(points)) {
    dataPath <- paste0(nc.path, "cmems/l4Chl/globalUpdated/", fishyear[i], "/", fishmonth[i], "/", 
                       fishyear[i], fishmonth[i], fishday[i], 
                       "_cmems_obs-oc_glo_bgc-plankton_my_l4-gapfree-multi-4km_P1D.nc") 
    
    # There are a few files missing for some earlier years (before ~ 2000)
    if(file.exists(dataPath)) {
      dataFile <- nc_open(dataPath) # 1-2 secs
    } else {
      next
    }
    # print(dataFile)
    lat <- ncvar_get(dataFile, "lat")
    lon <- ncvar_get(dataFile, "lon") # degrees west: -180 through 180
    lon360 <- ifelse(lon < 0, lon + 360, lon) # degrees east
    # tim <- ncvar_get(dataFile, 'time')
    nrows <- length(lon)
    ncols <- length(lat)
    
    #  Determine the pixels needed for fish location, if lat/long are not NA
    if (!is.na(points$lat[i]) & !is.na(points$lon360[i])) {
      # Calculate the location of the central pixel
      c <- which.min(abs(lon360 - points$lon360[i]))
      r <- which.min(abs(lat - points$lat[i]))
      # Use desired.diameter to calculate a pixel.diamater, and then a radius
      # If diameter is odd, apply radius around central pixel
      # If diameter is even, just trim off west/south to keep diameter correct
      pixel.diameter <- round(desired.diameter / 0.0416667, 0) 
      if((pixel.diameter %% 2) == 0) { # i.e. is an even number
        pixel.radius <- pixel.diameter / 2 
        c_low <- c - (pixel.radius - 1)
        c_up <- c + pixel.radius
        r_low <- r - (pixel.radius - 1)
        r_up <- r + pixel.radius
      } else { # i.e. is an odd number
        pixel.radius <- (pixel.diameter - 1) / 2 
        c_low <- c - pixel.radius
        c_up <- c + pixel.radius
        r_low <- r - pixel.radius
        r_up <- r + pixel.radius
      }
      
      # Calculate the number of rows and columns to extract
      numcols = abs(c_up - c_low) + 1 
      numrows = abs(r_up - r_low) + 1
      
      # Extract: Fix for if extraction box crosses the dateline (-180)
      if(c_low < 1) {
        data.var1 <- ncvar_get(dataFile, "CHL", start = c(1, r_low, 1), # just east of dateline
                              count = c(c_up, numrows, 1), verbose = FALSE) 
        data.var2 <- ncvar_get(dataFile, "CHL", start = c((8640 + c_low), r_low, 1), # just west of dateline
                               count = c(abs(c_low) + 1, numrows, 1), verbose = FALSE) 
        data.var <- rbind(data.var, data.var2)
      } else {
        data.var <- ncvar_get(dataFile, "CHL", start = c(c_low, r_low, 1), 
                              count = c(numcols, numrows, 1), verbose = FALSE)
      }
      
      # Add the correct column to points to receive values
      if(!is.na(mean(data.var, na.rm = TRUE))) {
        if (func == "mean") {
          points[i, ncol(points)] <- mean(data.var, na.rm = TRUE)
        } else if (func == "sd") {
          points[i, ncol(points)] <- sd(data.var, na.rm = TRUE)
        }
      }
    nc_close(dataFile)
    if((i - 100) %% 100 == 0) {
      print(paste0(y, " points complete"))
    }
   } 
  } # end i loop
  return(points)
}