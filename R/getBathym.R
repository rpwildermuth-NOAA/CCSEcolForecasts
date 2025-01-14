###########################################################################################################
# Extract bathymetry at points from ERDDAP. Can be slow for large datasets 
# Requires an input data frame with columns labeled "lon" and "lat"
# Note! In response to past inconsistencies, desired.resolution has been replaced with desired.diameter
# Desired.diameter = 0.7 (degrees) results in a 7x7 (= 49 pixel) box around the central location
# Contact Barbara.Muhling@noaa.gov
###########################################################################################################

require(ncdf4)
require(lubridate)
require(rerddap)
require(rerddapXtracto)

###########################################################################################################
# # Testing function
# points <- data.frame(expand.grid("lon" = seq(-132, -122, by = 2), "lat" = seq(32, 40, by = 2),
#                                  "date" = seq(as.Date("2009-06-15"), as.Date("2017-06-15"), by = "month")))
# desired.diameter <- 0.7 # Note! This means a 0.7x0.7 degree box
# func = "mean" # "mean" or "sd"
# testbathym <- getBathym(points = points, desired.diameter = desired.diameter, func = func)

###########################################################################################################

getBathym <- function(points, desired.diameter, func = "mean") {
  fishlat <- points$lat
  # See what longitude is called, and whether it's in degrees east or west
  # I often use "lon360" to specify longitude in degrees east: if that's the case, just grab that
  if("lon360" %in% colnames(points)) {
    fishlon360 <- points$lon360
  } else if("lon" %in% colnames(points)) { # Or maybe there's a column called "lon", and we don't know how it's measured
    fishlon360 <- points$lon
  }
  # Convert to degrees east if needed
  fishlon360 <- ifelse(fishlon360 < 0, fishlon360 + 360, fishlon360) 
  
  # Get dataset info from ERDDAP
  bathInfo <- rerddap::info('etopo360', url = 'http://coastwatch.pfeg.noaa.gov/erddap/')
  # Extract bathymetry at points
  bathym <- rxtracto(bathInfo, parameter = 'altitude', xcoord = fishlon360, ycoord = fishlat, 
                     xlen = desired.diameter, ylen = desired.diameter, progress_bar = TRUE)
  bathym <- data.frame(bathym)
  
  # Create cols
  points$out <- NA
  
  # Output values
  if(func == "mean") {
    points$out <- bathym$mean.altitude 
    colnames(points)[ncol(points)] <- paste0("bathym_mean", "_", desired.diameter)
  } else if (func == "sd") {
    points$out <- bathym$stdev.altitude
    colnames(points)[ncol(points)] <- paste0("bathym_sd", "_", desired.diameter)
  }
  points$count <- bathym$n
  colnames(points)[ncol(points)] <- paste0("bathym_count", "_", desired.diameter)
  # Return
  return(points)
}
