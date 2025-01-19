# Author: Allison Cluett
# Jan 2025
# From https://github.com/allisoncluett/cefi_westcoast_dst/blob/master/regrid.py
# 
# import xarray as xr
# import matplotlib.pyplot as plt
# import xesmf
# 
# def retrieve_static(path ='/Volumes/Cluett2024/MOM6/ocean_daily.static.nc'):
# 
#     NEP_static_file =  path
#     ds_static = xr.open_dataset(NEP_static_file)
# 
#     ds_static['geolon_c'] = xr.DataArray(data=ds_static['geolon_c'][1:, 1:], dims=('yq', 'xq'))
#     ds_static['geolat_c'] = xr.DataArray(data=ds_static['geolat_c'][1:, 1:], dims=('yq', 'xq'))
#     ds_static['geolon'] = xr.DataArray(data=ds_static['geolon'], dims=('yh', 'xh'))
#     ds_static['geolat'] = xr.DataArray(data=ds_static['geolat'], dims=('yh', 'xh'))
#     ds_static['lon'] = xr.DataArray(data=ds_static['geolon'])
#     ds_static['lat'] = xr.DataArray(data=ds_static['geolat'])
# 
#     return ds_static
# 
# ds_static = retrieve_static()
# 
# ## Read in MOM6 file
# mom6_file = '/Volumes/Cluett2024/MOM6/NEP10k_082024/ocean_daily/5yr/ocean_daily.19930101-19971231.tos.nc'
# ds = xr.open_dataset(mom6_file)
# 
# ds['tos'][1,:,:].plot()
# plt.show()
# 
# # Assign coordinates and regrid to MOM6
# roms_file = '/Volumes/Cluett2024/ROMS/wcnrt_daily/wcnrt_sst_daily_20110102_20240630.nc'
# roms = xr.open_dataset(roms_file)
# 
# roms['geolon'] = xr.DataArray(data=roms['lon_rho'], dims=('eta_rho', 'xi_rho'))
# roms['geolat'] = xr.DataArray(data=roms['lat_rho'], dims=('eta_rho', 'xi_rho'))
# 
# roms['sst'][1,:,:].plot()
# plt.show()
# 
# roms = roms.assign_coords({'lon': roms['geolon'],'lat': roms['geolat']})
# 
# roms_to_mom = xesmf.Regridder(roms, ds_static, method='bilinear', unmapped_to_nan=True)
# mom_to_roms = xesmf.Regridder(ds_static, roms, method='bilinear', unmapped_to_nan=True)
# 
# #mom_to_roms(ds['tos'][1,:,:]).plot()
# rg = mom_to_roms(ds['tos'])
# rg[1,:,:].plot()
# 
# plt.show()
# 
# #mom6_regrid = (mom_to_roms(ds)).to_netcdf(path='/Users/allisoncluett/Desktop')
