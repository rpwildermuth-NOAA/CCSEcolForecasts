#####################################################################################
# Extract MOM6 SST or chl
# myVar so far can be "chlos" (surface chlorophyll) or "tos" (SST)
# timestep can be "daily" or "monthly" (I don't have daily chl though, just SST)
# Note that unlike ROMS, MOM6 uses degrees east for longitude
# There's a bit of a kludgy fix where if you have any longitudes < 0, I assume they're degrees west
# and convert them. Input "points" need columns labelled "lat", "lon" and/or "lon360", and "date"
# There is some annoying fussing with filenames: this will likely not be needed in future once 
# MOM6 output filenames become more standardized!
# My netcdfs are in two subfolders of the "nc.path", "daily" and "monthly": code assumes yours are too
# Note I changed original filenames from what Allison sent me as they weren't clear
# I have uploaded them here: https://drive.google.com/drive/u/1/folders/1Jl00QqeG5okFOcxyGzhfaoCbWUj3Px38
# (Barb note to self: this is a mirror of script in xtractoLocal)
#####################################################################################

# # Test function
# points <- data.frame(expand.grid("lon" = seq(-134, -126, by = 2), "lat" = seq(30, 40, by = 2), 
#                                  "date" = seq(as.Date("1998-06-15"), as.Date("2017-06-15"), by = "year")))
# myVar <- "tos"
# pixel.radius <- 5
# timestep <- "daily"
# nc.path <- "F:/mom6/"
# test <- getMOM6(points = points, myVar = myVar, pixel.radius = pixel.radius, timestep = timestep, nc.path = nc.path)

getMOM6 <- function(points, myVar, pixel.radius, timestep, nc.path = "F:/mom6/") {
  library(lubridate)
  library(ncdf4)
  # Define biological observation locations and dates
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
  # Define date and year
  fishdate <- as.Date(points$date)
  fishyear <- year(fishdate)
  
  # Define output data
  out.data <- points
  
  # Will need the MOM6 grid files to get lon/lat locations (they aren't in the environmental netcdfs)
  grd <- nc_open(paste0(nc.path, "ocean_daily.static.nc"))
  geolon <- ncvar_get(grd, "geolon") 
  geolat <- ncvar_get(grd, "geolat")
  # Gives us a list of all possible lat/lons
  lonLong <- reshape2::melt(geolon, value.name = "lon")
  latLong <- reshape2::melt(geolat, value.name = "lat")
  ll <- data.frame(cbind("lon" = lonLong$lon, "lat" = latLong$lat))
  ll <- subset(ll, !is.na(lon) & !is.na(lat))
  nc_close(grd)
  
  # Remove cnk if testing code, causes issues!
  suppressWarnings(rm(cnk))
  
  ##########################################################################################################
  # Now extract environmental variable for all data points
  for(y in 1:nrow(points)) {
    # Skip if outside current MOM6 range (which may change in future!)
    if(fishyear[y] < 1998 | fishyear[y] > 2017 | is.na(fishlat[y] | is.na(fishlon360[y]))) {
      next
    }
    
    # First need to define year "chunk" name: MOM6 netcdfs currently provided as files containing 5 years each, part of filename
    # This also varies depending on if file contains daily or monthly outputs
    # If this is the first point to extract, we need to define the relevant netcdf and open it
    if(!exists("cnk")) {
      if(timestep == "daily") {
        cnk <- ifelse(fishyear[y] <= 2002, "19980101-20021231", ifelse(fishyear[y] <= 2007, "20030101-20071231", 
                                                            ifelse(fishyear[y] <= 2012, "20080101-20121231", "20130101-20171231")))
      } else if (timestep == "monthly") {
        cnk <- ifelse(fishyear[y] <= 2002, "199801-200212", ifelse(fishyear[y] <= 2007, "200301-200712", 
                                                                ifelse(fishyear[y] <= 2012, "200801-201212", "201301-201712")))
      }
      
      # Define the filename and open it
      dataPath <- paste0(nc.path, timestep, "/ocean_", timestep, ".", cnk, ".", myVar, ".nc")
      dataFile <- nc_open(dataPath) 
      
    # If this isn't the first point to extract, we can first check if we can just use the netcdf that's currently open
    } else {
      lastCnk <- cnk
      # Calculate the time chunk needed
      if(timestep == "daily") {
        cnk <- ifelse(fishyear[y] <= 2002, "19980101-20021231", ifelse(fishyear[y] <= 2007, "20030101-20071231", 
                                                                ifelse(fishyear[y] <= 2012, "20080101-20121231", "20130101-20171231")))
      } else if (timestep == "monthly") {
        cnk <- ifelse(fishyear[y] <= 2002, "199801-200212", ifelse(fishyear[y] <= 2007, "200301-200712", 
                                                                ifelse(fishyear[y] <= 2012, "200801-201212", "201301-201712")))
      }
      # Open new file, only if needed
      if(lastCnk != cnk) {
        nc_close(dataFile)
        dataPath <- paste0(nc.path, timestep, "/ocean_", timestep, ".", cnk, ".", myVar, ".nc")
        dataFile <- nc_open(dataPath) 
      }
    }
    tim <- ncvar_get(dataFile, 'time') ##### error
    timDate <- as.Date(as.POSIXct(tim * 86400, origin = "1993-01-01"))
    
    ##########################################################################################################
    # Grid is not square, so finding closest pixel by lon/lat separately (as we do for ROMS) doesn't work.
    # Instead, we need to minimize distances (in degrees) simultaneously
    cr <- which((abs(geolon - fishlon360[y]) + abs(geolat - fishlat[y])) == 
                  min((abs(geolon - fishlon360[y]) + abs(geolat - fishlat[y])), na.rm = T), arr.ind = TRUE)
    # Find date closest to desired one
    t <- which(abs(timDate - fishdate[y]) == min(abs(timDate - fishdate[y])))
    
    # Extract variable: first just at point 
    if(pixel.radius == 0) {
      envExtract <- ncvar_get(dataFile, myVar, start = c(cr[1], cr[2], t), count = c(1, 1, 1), verbose = FALSE) 
      envSD <- NA
    } else {
      # Or within a window defined by the pixel radius
      # Just using start and count as we usually do with ROMS won't work. 
      # For now, just defining points within pixel window and then looping through them, although is very inefficient!
      cr1 <- seq(from = cr[1] - pixel.radius, to = cr[1] + pixel.radius, by = 1)
      cr2 <- seq(from = cr[2] - pixel.radius, to = cr[2] + pixel.radius, by = 1)
      out <- data.frame(cbind("cr1" = cr1, "cr2" = cr2, "extract" = NA))
      for(i in 1:nrow(out)) {
        out$extract[i] <- ncvar_get(dataFile, myVar, start = c(out$cr1[i], out$cr2[i], t), count = c(1, 1, 1), verbose = FALSE)   
      }
      envExtract <- mean(out$extract, na.rm = TRUE)
      envSD <- sd(out$extract, na.rm = TRUE)
    }
    
    ##########################################################################################################
    # Save to out.data. If y is 1 (i.e. first row of out.data), must make new column before filling it
    if (y == 1) {
      out.data$var <- NA
      out.data$varsd <- NA
    }
    out.data$var[y] <- envExtract
    out.data$varsd[y] <- envSD 
    
    if((y - 100) %% 100 == 0) {
      print(paste0(y, " points complete"))
    }
  } 
  nc_close(dataFile)
  
  # Colnames use myVar and timestep 
  colnames(out.data)[ncol(out.data) - 1] <- paste0(myVar, "_", timestep)
  colnames(out.data)[ncol(out.data)] <- paste0(myVar, "_", timestep, "_sd")
  return(out.data)
}
