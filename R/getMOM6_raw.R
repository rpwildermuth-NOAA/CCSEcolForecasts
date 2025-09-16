###########################################################################################################
# Extract MOM6 from the native grid
# timestep is daily
# Note I changed original filenames from what Allison sent me as they weren't clear
# I obtained the nc files from here: https://psl.noaa.gov/cefi_portal/
# # Script created by Barbara Muhling and modified by Nerea Lezama-Ochoa
###########################################################################################################

#libraries
require(lubridate)
require(ncdf4)
require(dplyr)



###########################################################################################################
getMOM6_raw <- function(points, varName, desired.diameter, timestep, func = "mean", nc.path = "C:\\Users\\nereo\\Dropbox (Personal)\\Nerea\\NOAA\\PROJECTS & COLLABORATIONS\\PROJECTS\\Forecast\\Swordfish_forecast\\1_extract_variables\\MOM6/") {
  # Define biological observation locations and dates
  fishlat <- swor$lat
  # See what longitude is called, and whether it's in degrees east or west
  # I often use "lon360" to specify longitude in degrees east: if that's the case, just grab that
  if("lon360" %in% colnames(swor)) {
    fishlon360 <- swor$lon360
  } else if("lon" %in% colnames(swor)) { # Or maybe there's a column called "lon", and we don't know how it's measured
    fishlon360 <- swor$lon
  }
  # Convert to degrees east if needed
  fishlon360 <- ifelse(fishlon360 < 0, fishlon360 + 360, fishlon360) 
  # Define date and year
  fishdate <- as.Date(swor$date)
  fishyear <- year(fishdate)
  
  # Define output data
  out.data <- swor
  
  # Define output colnames
  if(func == "mean") {
    out.data$mn <- NA
    colnames(out.data)[ncol(out.data)] <- paste0(varName, "_", "mean", "_", desired.diameter, "_", timestep)
    out.data$count <- NA
    colnames(out.data)[ncol(out.data)] <- paste0(varName, "_", "count", "_", desired.diameter, "_", timestep)
  }
  if(func == "sd") {
    out.data$sd <- NA
    colnames(out.data)[ncol(out.data)] <- paste0(varName, "_", "sd", "_", desired.diameter, "_", timestep)
    out.data$count <- NA
    colnames(out.data)[ncol(out.data)] <- paste0(varName, "_", "count", "_", desired.diameter, "_", timestep)
  }
  
  # Will need the MOM6 grid files to get lon/lat locations (they aren't in the environmental netcdfs)
  grd <- nc_open(paste0(nc.path, "ocean_daily.static.nc"))
  geolon <- ncvar_get(grd, "geolon") 
  geolat <- ncvar_get(grd, "geolat")
  nc_close(grd)
  
  # Remove cnk if testing code, causes issues
  suppressWarnings(rm(cnk))
  
  ##########################################################################################################
  # Now extract environmental variable for all data points
  for(y in 1:nrow(points)) {
    # Skip if outside current MOM6 range (which may change in future!)
    if(fishyear[y] < 1993 | fishyear[y] > 2017 | is.na(fishlat[y] | is.na(fishlon360[y]))) {
      next
    }
    
    # First need to define year "chunk" name: MOM6 netcdfs currently provided as files containing 5 years each, part of filename
    # This also varies depending on if file contains daily or monthly outputs
    # If this is the first point to extract, we need to define the relevant netcdf and open it
    if(!exists("cnk")) {
      if(timestep == "daily") {
        cnk <- ifelse(fishyear[y] <= 2019, "19930101-20191231", ifelse(fishyear[y] <= 2024, "20200101-20241231", 
                                                                       ifelse(fishyear[y] <= 2029, "20250101-20291231", "20250101-20291231")))
      } 
      
      # Define the filename and open it
      dataPath <- paste0(nc.path, timestep, "/ocean_", timestep, ".", cnk, ".", varName, ".nc")
      dataFile <- nc_open(dataPath) 
      
      # If this isn't the first point to extract, we can first check if we can just use the netcdf that's currently open
    } else {
      lastCnk <- cnk
      # Calculate the time chunk needed
      if(timestep == "daily") {
        cnk <- ifelse(fishyear[y] <= 2019, "19930101-20191231", ifelse(fishyear[y] <= 2024, "20200101-20241231", 
                                                                       ifelse(fishyear[y] <= 2029, "20250101-20291231", "20250101-20291231")))
      } 
      # Open new file, only if needed
      if(lastCnk != cnk) {
        nc_close(dataFile)
        dataPath <- paste0(nc.path, timestep, "/ocean_", timestep, ".", cnk, ".", varName, ".nc")
        dataFile <- nc_open(dataPath) 
      }
    }
    tim <- ncvar_get(dataFile, 'time') 
    timDate <- as.Date(as.POSIXct(tim * 86400, origin = "1993-01-01"))
    # Find date closest to desired one
    t <- which(abs(timDate - fishdate[y]) == min(abs(timDate - fishdate[y])))
    
    ##########################################################################################################
    # MOM6 grid is not square, so finding closest pixel by lon/lat separately (as we do for ROMS) doesn't work.
    # Instead, we can minimize distances (in degrees) simultaneously
    # First check to make sure point isn't far outside MOM6 domain
    minDist <- min((abs(geolon - fishlon360[y]) + abs(geolat - fishlat[y])), na.rm = T)
    if(minDist > 0.2) {
      next
    }
    # If not, calculate closest grid point
    cr <- which((abs(geolon - fishlon360[y]) + abs(geolat - fishlat[y])) == 
                  min((abs(geolon - fishlon360[y]) + abs(geolat - fishlat[y])), na.rm = T), arr.ind = TRUE)
    
    # Defining a pixel extraction window is tricky: the number of pixels corresponding to the window (in degrees)
    # will vary with location
    # First define the min and max of lon and lat in the desired pixel extraction box
    lonmin <- fishlon360[y] - (desired.diameter / 2)
    lonmax <- fishlon360[y] + (desired.diameter / 2)
    latmin <- fishlat[y] - (desired.diameter / 2)
    latmax <- fishlat[y] + (desired.diameter / 2)
    
    # If range only covers one pixel, extract at that pixel
    if(lonmin == lonmax & latmin == latmax) { 
      envExtract <- ncvar_get(dataFile, varName, start = c(cr[1], cr[2], t), count = c(1, 1, 1), verbose = FALSE)
      envPixels <- 1
    } else {
      # Or within a window defined by the desired.diameter (degrees) 
      # Just using start and count as we usually do with ROMS won't work
      # MOM6 cell sizes also change substantially with longitude and latitude
      # So for each point, we need to find the closest points defining the desired.diameter, and then extract those values
      # First define all lon and lat points within the desired range
      lonRange <- data.frame(which(geolon <= lonmax & geolon >= lonmin, arr.ind = TRUE))
      latRange <- data.frame(which(geolat <= latmax & geolat >= latmin, arr.ind = TRUE))
      # Inner join to get only row/cols in geolon/geolat that meet BOTH conditions
      boxRange <- inner_join(lonRange, latRange, by = c("row", "col"))
      # If the desired.diameter is small, it could only cover 1 pixel anyway, so extract at that point
      if(nrow(boxRange) <= 1) {
        envExtract <- ncvar_get(dataFile, varName, start = c(cr[1], cr[2], t), count = c(1, 1, 1), verbose = FALSE)
        envPixels <- 1
      } else {
        # plot(boxRange$col, boxRange$row) # Can see pixel extraction box is not square
        # Fastest approach may be to extract all pixels in square covering boxRange, and then drop those not needed
        lonStart <- min(boxRange$row)
        lonCount <- (max(boxRange$row) - lonStart) + 1
        latStart <- min(boxRange$col)
        latCount <- (max(boxRange$col) - latStart) + 1
        envExtractBox <- ncvar_get(dataFile, varName, start = c(lonStart, latStart, t), count = c(lonCount, latCount, 1), 
                                   verbose = FALSE, collapse_degen = FALSE) # Last arg stops singleton/degen dimensions being dropped 
        # Drop extra time dimension
        envExtractBox <- as.matrix(envExtractBox[, , 1], drop = FALSE, byrow = TRUE)
        # Special case if envExtractBox has 2 rows and 1 col: as.matrix automatically transposes it, which we don't want
        if(length(unique(boxRange$row)) == 1 & length(unique(boxRange$col)) == 2) {
          envExtractBox <- t(envExtractBox)
        }
        # Subset this rectangular matrix using boxRange indices
        dimnames(envExtractBox) <- list(row = seq(min(boxRange$row), max(boxRange$row)),
                                        col = seq(min(boxRange$col), max(boxRange$col))) 
        envExtractLong <- as.data.frame.table(envExtractBox, stringsAsFactors = FALSE)
        envExtractLong$row <- as.numeric(envExtractLong$row)
        envExtractLong$col <- as.numeric(envExtractLong$col)
        colnames(envExtractLong)[3] <- varName
        envExtract <- inner_join(envExtractLong, boxRange, by = c("row", "col"))
        # Also record number of non-NA pixels: this will change with latitude
        envPixels <- sum(!is.na(envExtract[, 3])) 
      }
    }
    
    ##########################################################################################################
    # Save values to out.data
    if (func == "mean") {
      if(nrow(envExtract) == 1) { # single point
        out.data[y, ncol(out.data) - 1] <- envExtract
      } else { # mean of matrix
        out.data[y, ncol(out.data) - 1] <- mean(envExtract[, 3], na.rm = TRUE)
      }
      out.data[y, ncol(out.data)] <- envPixels 
    } else if (func == "sd") {
      if(nrow(envExtract) == 1) { # single point
        out.data[y, ncol(out.data) - 1] <- NA
      } else { # sd of matrix
        out.data[y, ncol(out.data) - 1] <- sd(envExtract[, 3], na.rm = TRUE)
      }
      out.data[y, ncol(out.data)] <- envPixels
    }
    
    if((y - 100) %% 100 == 0) {
      print(paste0(y, " points complete"))
    }
  } 
  nc_close(dataFile)
  return(out.data)
}
