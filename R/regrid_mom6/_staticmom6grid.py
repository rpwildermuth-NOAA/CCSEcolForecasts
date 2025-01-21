# Python script for renaming lat/lon variables in `ocean_daily.static.nc`
# Code authors: Allison Cluett/Jessica Bolin
# January 2025

#### DO NOT RUN THIS SCRIPT! #####
# Jessie has already ran this script in PyCharm, which produced ds_static.nc. 
# This is the empty MOM6 grid we use for regridding. I have just saved this script 
# for posterity/reproducibility. 

## Import modules (libraries)
#import xarray as xr #for opening netcdfs

# Open static ocean model grid
#path = "data/ocean_daily.static.nc"
#NEP_static_file = path
#ds_static = xr.open_dataset(NEP_static_file)

# Rename lat/lon variables. Ask Allison for help understanding python slicing conventions 
#ds_static['geolon_c'] = xr.DataArray(data=ds_static['geolon_c'][1:, 1:], dims=('yq', 'xq'))
#ds_static['geolat_c'] = xr.DataArray(data=ds_static['geolat_c'][1:, 1:], dims=('yq', 'xq'))
#ds_static['geolon'] = xr.DataArray(data=ds_static['geolon'], dims=('yh', 'xh'))
#ds_static['geolat'] = xr.DataArray(data=ds_static['geolat'], dims=('yh', 'xh'))
#ds_static['lon'] = xr.DataArray(data=ds_static['geolon'])
#ds_static['lat'] = xr.DataArray(data=ds_static['geolat'])

# Write to data directory 
#ds_static.to_netcdf(path='data/data/ds_static.nc')
