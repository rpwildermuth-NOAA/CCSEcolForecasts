###############################################################################################################
# Calculate distance to nearest land in the eastern North Pacific
# Requires a coast shapefile, which is clunky. Can be very slow for a lot of points
# Points must have cols lon and lat, lon is in degW (i.e. does not work across dateline)
###############################################################################################################

require(sf)
require(nngeo) 

getDistLand <- function(points, coast) { 
  # Define projections
  epsg.2062 <- "+proj=lcc +lat_1=40 +lat_0=40 +lon_0=0 +k_0=0.9988085293 +x_0=600000 +y_0=600000 +a=6378298.3 +b=6356657.142669561 +pm=madrid +units=m +no_defs"
  wgs.84    <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  # st_crs(coast) # WGS 84
  coast.proj <- st_transform(coast, epsg.2062) 
  
  # Extract distance from nearest land. Can be slow
  toCalc <- points[c("lon", "lat")]
  toCalc <- st_as_sf(toCalc, coords = c("lon", "lat"), crs = wgs.84)
  toCalcProj <- st_transform(toCalc, epsg.2062)
  dists <- st_nn(toCalcProj, coast.proj, k = 1, parallel = 5, returnDist = TRUE)
  dists.vec <- unlist(dists$dist)
  points$distLand <- dists.vec
  return(points)
}