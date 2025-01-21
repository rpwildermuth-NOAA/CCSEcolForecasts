# Function for regridding MOM6 using R and Python 
# Authors: Jessica Bolin and Allison Cluett
# January 2025

# You will need to add the GDrive folder `regrid_example_jb` to the /data 
# directory in the Github repo.
# I've added data/regrid_example_jb to the .gitignore. 

# Dependencies ------------------------------------------------------------

library(reticulate) 
library(terra)
library(maps)

packageVersion("reticulate") #‘1.40.0’
packageVersion("terra") #‘1.7.78’
packageVersion("maps") #‘3.4.2’

pth <- "/Users/admin/Documents/GitHub/CCSEcolForecasts/data/regrid_example_jb/"

# Set up python environment -----------------------------------------------

use_python("/opt/anaconda3/envs/xesmf_env_2/bin/python")
sys <- import("sys")            # for checking python version/env 
xarray <- import("xarray")      # for opening netcdfs
matplotlib <- import("matplotlib")  # for visualization
plt <- matplotlib$pyplot        # for visualization using matplotlib
xesmf <- import("xesmf")        # for regridding
sys$version # Jessie: "3.12.8 | packaged by conda-forge | (main, Dec  5 2024, 14:23:40) [Clang 18.1.8 ]"
sys$executable #Jessie: "/opt/anaconda3/envs/xesmf_env/bin/python"


# Files -------------------------------------------------------------------

NEP_static_file <- paste0(pth, "ds_static.nc") # Read in static MOM6 grid 
ds_static <- xarray$open_dataset(NEP_static_file)

roms_file <- paste0(pth, "gfdltv_1990_temp.nc") # ROMS grid 
roms <- xarray$open_dataset(roms_file)
roms <- roms$assign_coords(list(lon = roms$longitude, lat = roms$latitude)) 

fileys <- list.files(pth)[grep("ocean_monthly", list.files(pth))]  # MOMS files 
fileys <- fileys[!grepl("_RG.nc", fileys)] #exclude any regridded files that already exist in directory
fileys

# Function ----------------------------------------------------------------

regrid_MOM6 <- function(mom6_file) {

  filename <- gsub(pth, "", mom6_file) 
  filename <- gsub(".nc", "_RG.nc", filename) 
  
  if (!file.exists(paste0(pth, filename))) {
    
  ds <- xarray$open_dataset(mom6_file)
  mom_to_roms <- xesmf$Regridder(ds_static, roms, method = 'bilinear', unmapped_to_nan = TRUE)
  mom6_regrid <- (mom_to_roms(ds))$to_netcdf(path = paste0(pth, filename))
  } else { 
    print(paste0("skipping iteration ", fileys[i])) 
    
    } # if then
} # function

# Run function ------------------------------------------------------------

for (i in 1:length(fileys)) {
  regrid_MOM6(mom6_file = paste0(pth, fileys[i]))
}


# Visusalise, check it worked using R ---------------------------------------------

rr <- rast(paste0(pth, "ocean_monthly.201301-201712.tos_RG.nc"))
r <- rr[[50]] 
plot(r); maps::map("world", add = T) # *muscle emoji* 
r #3km resolution, WGGS84 projection, same as ROMS, looks promising

rr <- rast(paste0(pth, "ocean_monthly.199801-200212.chlos_RG.nc"))
r <- rr[[50]] 
plot(r); maps::map("world", add = T) # *muscle emoji* 
