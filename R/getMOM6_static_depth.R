getMOM6_static_depth <- function(points, varName = "deptho", desired.diameter,func = "mean", nc.file = "C:\\Users\\nereo\\Dropbox (Personal)\\Nerea\\NOAA\\PROJECTS & COLLABORATIONS\\PROJECTS\\Forecast\\Swordfish_forecast\\1_extract_variables\\MOM6\\daily/ocean_daily.19930101-20191231.deptho.nc") {
  require(ncdf4)
  require(dplyr)
  
  points=swor
  # Convert lon to 0–360 if needed
  if ("lon360" %in% colnames(points)) {
    points$lon360 <- points$lon360
  } else if ("lon" %in% colnames(points)) {
    points$lon360 <- ifelse(points$lon < 0, points$lon + 360, points$lon)
  }
  
  # Open the static NetCDF file and get the grid
  nc <- nc_open(nc.file)
  geolon <- ncvar_get(nc, "geolon")
  geolat <- ncvar_get(nc, "geolat")
  depth_grid <- ncvar_get(nc, varName)
  nc_close(nc)
  
  # Initialize output
  points[[paste0(varName, "_value")]] <- NA
  
  # Extract depth at nearest grid cell
  for (i in 1:nrow(points)) {
    lon <- points$lon360[i]
    lat <- points$lat[i]
    
    if (is.na(lon) | is.na(lat)) next
    
    # Find nearest point in grid
    dist <- abs(geolon - lon) + abs(geolat - lat)
    if (min(dist, na.rm = TRUE) > 0.2) next  # skip points outside domain
    
    nearest <- which(dist == min(dist, na.rm = TRUE), arr.ind = TRUE)
    points[[paste0(varName, "_value")]][i] <- depth_grid[nearest[1], nearest[2]]
  }
  
  return(points)
}
