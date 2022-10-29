#Step 1: Start with import statements (this tells python which packages 
#you will actually use in your code and names 
#them, so we will call xarray : xr)
import sys
import xarray as xr
import pandas as pd
import numpy as np
import glob
import collections

#%%
##Step 2: Read in our data from a netcdf file of rainfall over Ethiopia
#(edit the path to the file correctly depending on where you saved it) 
#and select the variable ‘precip’. Now we have a data array called da and
#we can print it to see what we have:
# da = xr.open_dataarray("/Volumes/USBkey32GB/water_ex/chirps-v2.0.monthly.nc")
# da = da.rename({'latitude':'lat','longitude':'lon'})
# da = da.sel(lat=slice(3,15),lon=slice(25,45))
# da.to_netcdf('/Volumes/USBkey32GB/water_ex/reg_chirps_monthly.nc')

da = xr.open_dataset("/Volumes/USBkey64GB/health_ex/tamsat_daily_ethiopia.nc")['rfe_filled']
#convert from mm/month to mm/day
da = da.sel(time=slice('1981-01-01','2020-12-31'))
# print data array
print(da)

tda = xr.open_mfdataset("/Volumes/USBkey64GB/water_ex/CRU_TS4_06/cru_ts4.06.*.tmp.dat.nc")
print(tda)
tda = tda['tmp']

eda = xr.open_dataset("/Volumes/USBkey64GB/health_ex/ERA5_surface.nc")
print(eda)




#%%
# rasterarray.rio.write_crs("epsg:4326", inplace=True)

# data = gpd.read_file('../Awash/Awash_basin_border.shp')
# print(data)
# print(data.keys())
# print(data['OBJECTID'])
# print(data['BASIN_ID'])
# print(data['geometry'])
# print("awash crs", data.crs)
# new_seas_anom = seas_anom.rio.clip(data.geometry.apply(mapping),data.crs)
#%%

