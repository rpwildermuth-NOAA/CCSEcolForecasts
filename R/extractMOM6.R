# code to extract MOM6 variables for CPS SDM models

library("ncdf4")

# depth and lat/long
# Specify the OPeNDAP server URL (using regular grid output)
url_static <-"http://psl.noaa.gov/thredds/dodsC/Projects/CEFI/regional_mom6/cefi_portal/northeast_pacific/full_domain/hindcast/daily/raw/r20250912/ocean_static.nc"
ncopendap_static <- nc_open(url_static)
depthNEP <- ncvar_get(ncopendap_static, "deptho")
lon <- ncvar_get(ncopendap_static, "None")
lat <- ncvar_get(ncopendap_static, "None")

# SSH
url_ssh <-"http://psl.noaa.gov/thredds/dodsC/Projects/CEFI/regional_mom6/cefi_portal/northeast_pacific/full_domain/hindcast/daily/raw/r20250912/ssh.nep.full.hcast.daily.raw.r20250912.199301-202506.nc"
ncopendap <- nc_open(url_ssh)
sshNEP <- ncvar_get(ncopendap, "ssh")

# SST
url_sst <-"http://psl.noaa.gov/thredds/dodsC/Projects/CEFI/regional_mom6/cefi_portal/northeast_pacific/full_domain/hindcast/daily/raw/r20250912/tos.nep.full.hcast.daily.raw.r20250912.199301-202506.nc"
ncopendap <- nc_open(url_sst)
sstNEPdaily <- ncvar_get(ncopendap, "tos")
url_sst <-"http://psl.noaa.gov/thredds/dodsC/Projects/CEFI/regional_mom6/cefi_portal/northeast_pacific/full_domain/hindcast/monthly/raw/r20250912/tos.nep.full.hcast.monthly.raw.r20250912.199301-202506.nc"
ncopendap <- nc_open(url_sst)
sstNEPmonthly <- ncvar_get(ncopendap, "tos")
