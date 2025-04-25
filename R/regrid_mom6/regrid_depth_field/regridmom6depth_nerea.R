# Regrid MOM6 depth to 3km ROMS for Nerea
# Jessica Bolin
# April 2025
# Uses workflow here: https://jessicabolin.quarto.pub/how-to-regrid-mom6-using-r-and-python/

# Dependencies ------------------------------------------------------------

library(reticulate) # for using python in R
library(terra) #raster data
library(maps) # country outline

pth <- "/Users/admin/Documents/GitHub/CCSEcolForecasts/data/regrid_example_jb/"
use_python("/opt/anaconda3/envs/xesmf_env_2/bin/python")

sys <- import("sys")            # for checking python version/env 
xarray <- import("xarray")      # for opening netcdfs
matplotlib <- import("matplotlib")  # for visualization
plt <- matplotlib$pyplot        # for visualization using matplotlib
xesmf <- import("xesmf")        # for regridding (takes a while)

sys$version  #"3.12.8  (not 3.13 - yay)
sys$executable #[1] "/opt/anaconda3/envs/xesmf_env_2/bin/python"


# Depth -------------------------------------------------------------------

# Create MOM6 depth field that we want to regrid. Subset the depth field,
# and write it as a netcdf to working directory 
NEP_static_file <- paste0(pth, "ds_static.nc")
ds_static <- xarray$open_dataset(NEP_static_file)
depth <- ds_static$deptho
depth_ds <- depth$to_dataset(name = "deptho")
depth_ds$to_netcdf(paste0(pth, "depth/deptho_field.nc")) #nice

# Read in example ROMS/rast  ----------------------------------------------

# Read in a pre-determined raster wiht the extent, CRS, projection etc that I 
# want. Here I'm using Mer's 3km ROMS regridded to off California. 
#roms_file <- paste0("/Users/admin/Documents/GitHub/CCSEcolForecasts/data/regrid_example_jb/gfdltv_1990_temp.nc")
#roms <- xarray$open_dataset(roms_file) # read in with python

# UPDATE. Make an empty raster specified for Nerea
# resolution; 0.2883739, 0.08588319
# crs: +proj=longlat +datum=wgs84 +no_defs
# extent: 156.7806, 255.1161, 10.76609, 80.76089

romsy <- rast(res = c(0.2883739, 0.08588319),
       ext(156.7806, 255.1161, 10.76609, 80.76089))
crs(romsy) <-  "+proj=longlat +datum=WGS84 +no_defs"
plot(romsy, axes = T); maps::map("world", add= T, wrap = c(0, 360))
#terra::writeCDF(romsy, paste0(pth, "depth/empty_rast_nerea.nc"))

roms_file <- "empty_rast_nerea.nc"
roms <- xarray$open_dataset(paste0(pth, "depth/", roms_file)) # read in with python

# Read in depth file created previously -----------------------------------

# Read it in and add _RG to the filename for when we regrid and save to working
#directory
mom6_file = paste0(pth, "depth/deptho_field.nc")
filename <- gsub(".nc", "_RG.nc", mom6_file) 
filename #[1] "/Users/admin/Documents/GitHub/CCSEcolForecasts/data/regrid_example_jb/depth/deptho_field_RG.nc"

#Open as xarray dataset
ds <- xarray$open_dataset(mom6_file)

# Function to regrid ------------------------------------------------------

#Here we are using bilinear interpolation for regridding the MOM6 grid to the 3 km ROMS grid.
mom_to_roms <- xesmf$Regridder(ds_static, # MOM6 empty grid
                               roms, # ROMS grid
                               method = 'bilinear', 
                               unmapped_to_nan = TRUE)

# Run function and save netcdf to wroking directory
mom_to_roms(ds)$to_netcdf(path = paste0(filename))
list.files(paste0(pth, "depth")) #[1] "deptho_field_RG.nc"  "deptho_field.nc"     "empty_rast_nerea.nc"

# Check with terra --------------------------------------------------------

rast(filename) %>% plot(main = "MOM6 depth (m)")
maps::map("world", add = T, wrap = c(0, 360), fill = T, col = "grey") # NICE! I love when code works 
rast(filename)
# class       : SpatRaster 
# dimensions  : 815, 341, 1  (nrow, ncol, nlyr)
# resolution  : 0.2883739, 0.08588319  (x, y)
# extent      : 156.7806, 255.1161, 10.76609, 80.76089  (xmin, xmax, ymin, ymax)
# coord. ref. : lon/lat WGS 84 (CRS84) (OGC:CRS84) 
# source      : deptho_field_RG.nc 
# varname     : deptho 
# name        : deptho 

