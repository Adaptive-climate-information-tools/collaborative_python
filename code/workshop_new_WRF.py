#Start with import statements (this tells python which packages 
#you will actually use in your code and names 
#them, so we will call xarray : xr)
#import sys
import xarray as xr
import xesmf as xe
import matplotlib.pyplot as plt
#%%
#Regrid WRF data to a regular lat lon grid using XESMF
#This will make it much easier to compare to other datasets 
#and apply evaluation code to
#XESMF can easily deal with the Lambert Conformal grid used by WRF
#https://xesmf.readthedocs.io/en/latest/why.html?highlight=wrf#for-emerging-new-grid-types
path = '/Users/ellendyer/Library/Mobile Documents/com~apple~CloudDocs/1SHARED_WORK/Work/REACH/Workshop_conda_python/'
path_wrf = '/Volumes/passport_ellen/'
ds = xr.open_dataset(path_wrf+'wrfout_d02_2022102006',engine='netcdf4')
#Select variable that you want from WRF file
#to select one variable:
# ds = ds['RAINNC']
#to select multiple variables:
ds = ds[['RAINNC','T2']]
ds = ds.rename({'Time':'time','south_north':'lat', 'west_east':'lon'})
#sort out time before regridding - rename and sort which is required
#for most pandas and xarray time operations
ds = ds.assign_coords(time=ds.XTIME.values)
ds = ds.drop('XTIME')
ds = ds.sortby('time')

# --test plot to compare later---
#**if selected on variable
# ds.mean('time').plot(robust=True)
# plt.show()
#**if selected more than one variable can plot 
#**test for each
ds['T2'].mean('time').plot(vmin=280,vmax=300)
plt.show()
ds['RAINNC'].mean('time').plot(vmin=0,vmax=5)
plt.show()

#regridding based on this example:
#https://github.com/zxdawn/GEOSChem-python-tutorial/blob/master/Chapter04_regridding_WRFChem.ipynb
#Pick your new grid resolution and get lat and lon bounds from wrf dataset
resolution = 0.1 #can chose anything the same degree size as WRF grid or bigger
dst = ds.isel(time=0)
lon_min = dst.XLONG.min()
lon_max = dst.XLONG.max()
lat_min = dst.XLAT.min()
lat_max = dst.XLAT.max()
#new grid to interpolate WRF to
newgrid = xe.util.grid_2d(lon_min-resolution, lon_max+resolution, resolution,
                        lat_min-resolution, lat_max+resolution, resolution)
regridder = xe.Regridder(dst, newgrid, method='bilinear')
#Apply regridder
wrfrg = regridder(ds)

#Rename dimensions and organise coordinates to be easy 
#to do time operations with xarray and to be 
#easily compared with other gridded datasets
lats, lons = wrfrg.lat[:,0].values, wrfrg.lon[0,:].values
wrfrg = wrfrg.drop(['lat','lon'])
wrfrg = wrfrg.rename({'y':'lat', 'x':'lon'})
wrfrg = wrfrg.assign_coords(lat=lats,lon=lons)

# --test plot to compare ---
#**if selected on variable
# wrfrg.mean('time').plot(robust=True)
# plt.show()
#**if selected more than one variable can plot 
#**test for each
wrfrg['T2'].mean('time').plot(vmin=280,vmax=300)
plt.show()
wrfrg['RAINNC'].mean('time').plot(vmin=0,vmax=5)
plt.show()

#Resample to daily data using mean - make sure this is the right
#operation given that you are working with accumulation output!
wrf_day = wrfrg.resample(time='D').mean()

#Save daily data to netcdf file
wrf_day.to_netcdf(path+'files/wrf_daily_out.nc')

