#Step 1: Start with import statements (this tells python which packages 
#you will actually use in your code and names 
#them, so we will call xarray : xr)

import xarray as xr
import numpy as np
#Import for Step 3
from matplotlib import pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature
#Import for Step 4
from scipy.stats import spearmanr, pearsonr
#Import for Step 6
import geopandas as gpd
from shapely.geometry import mapping
import rioxarray
import sys

#%%
##Step 2: Read in our data from a netcdf file of rainfall over Ethiopia
#(edit the path to the file correctly depending on where you saved it) 
#and select the variable ‘precip’. Now we have a data array called da and
#we can print it to see what we have:
path = '/Users/ellendyer/Library/Mobile Documents/com~apple~CloudDocs/1SHARED_WORK/Work/REACH/Workshop_conda_python/'
da = xr.open_dataarray(path+'workshop_setup/workshop_chirps.nc')
#convert from mm/month to mm/day
da = da/da.time.dt.daysinmonth
da = da.sel(time=slice('1981-01-01','2019-12-31'))
# print data array
print(da)
hsst = xr.open_dataset(path+'files/HadISST_sst.nc')['sst']
#Read in datafile if it is in the same directory
#hsst = xr.open_dataset('HadISST_sst.nc')['sst']
hsst = hsst.sel(time=slice('1981-01-01','2019-12-31'))
print(hsst)
# convert degree C to Kelvin
hsst = hsst + 273.15

#%%
##Step 3: Select a region and time slice from ts timeseries (you can use date strings
#with xarray). The .sel slice method chooses the nearest point so
#you don't need to know exact bounds

# Ethiopian rainfall
p_eth = da.sel(lat=slice(2,22),lon=slice(30,52))

# Global SSTs
gsst = hsst.sel(latitude=slice(60,-60),longitude=slice(-180,180))

# Western Indian Ocean SSTs
wsst = gsst.sel(latitude=slice(25,0),longitude=slice(45,70))
wsst_mean = wsst.mean(dim=('latitude','longitude'))

# Eastern Indian Ocean SSTs
esst = gsst.sel(latitude=slice(10,-10),longitude=slice(90,110))
esst_mean = esst.mean(dim=('latitude','longitude'))

# Nino 3.4 SSTs
nsst = gsst.sel(latitude=slice(5,-5),longitude=slice(-170,-120))
nsst_mean = nsst.mean(dim=('latitude','longitude'))

# Indian Ocean Dipole Index
iodw = gsst.sel(latitude=slice(10,-10),longitude=slice(50,70)).mean(dim=('latitude','longitude'))
iode = gsst.sel(latitude=slice(0,-10),longitude=slice(90,110)).mean(dim=('latitude','longitude'))
iodi_mean = iodw - iode

#YOUR INDEX HERE....

#%%
##Step 4: Calculate the seasonal mean of your index and select composite years 
#You can choose the number of years and the characteristic of these years

#Dipole mode index (Indian Ocean Dipole) example:
sst_index = iodi_mean
months = [3,4,5]
seas = 'MAM'
sst_seas_index = (sst_index.sel(time=sst_index.time.dt.month.isin(months))).groupby('time.year').mean('time')
#This sorts the array of dipole mode indices - it will go from smallest value (negative) to largest value (positive)
sort_index = sst_seas_index.sortby(sst_seas_index)
#Select the years for the composites (I am choosing the strongest positive and negative dipole years)
neg_comp_y = sort_index[0:6]['year']
pos_comp_y = sort_index[-6:]['year']

#%%
##Step 5: Calculate a seasonal timeseries of the variable you are interested in compositing
#In this example the variable is rainfall
#Create a seasonal average timeseries and a long-term seasonal average
#months = [7,8,9]
#seas = 'JAS'
rain_seas = da.sel(time=da.time.dt.month.isin(months)).groupby('time.year').mean('time')
rain_seas_mean = rain_seas.mean('year')
#Select the negative and positive composite years from the seasonal timeseries and create
#a composite mean
rain_comp_neg = rain_seas.sel(year=rain_seas.year.isin(neg_comp_y)).mean('year')
rain_comp_pos = rain_seas.sel(year=rain_seas.year.isin(pos_comp_y)).mean('year')
#Subtract the long term mean from the composite mean
#to get a composite anomaly
comp_neg_anom = rain_comp_neg - rain_seas_mean
comp_pos_anom = rain_comp_pos - rain_seas_mean

#%%
##Step 6: Plot the negative composite anomaly for all of Ethiopia
#Using the same plotting setup that we used for correlations 
#and station/shapefile code
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()
ax.add_feature(cartopy.feature.BORDERS,lw=2)
ax.add_feature(cartopy.feature.RIVERS,lw=2)
comp_neg_anom.plot(ax=ax,cmap = plt.cm.BrBG,transform=ccrs.PlateCarree()
               #,vmin=-0.8,vmax=0.8
               ,extend='both'
               ,robust=True,cbar_kwargs={'label': "mm/day"})
plt.title(seas+' negative IOD composite rainfall')
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='black', alpha=0.5, linestyle='dotted')
# for extent the order is  [West,East,South,North]
ax.set_extent([31, 48, 3, 15])
gl.top_labels = False
gl.right_labels = False
plt.savefig(path+'plots/C_1.png',bbox_inches='tight',dpi=200)
#Save in the same directory as code
#plt.savefig('C_1.png',bbox_inches='tight',dpi=200)
plt.show()
plt.clf()

#%%
##Step 7: Plot the positive composite anomaly for the Awash basin

#Make sure the composite anomaly has the specific coordinate 
#system of the shapefile
comp_pos_anom.rio.write_crs("epsg:4326", inplace=True)
#Read in Awash basin shapefile
data = gpd.read_file(path+'Awash/Awash_basin_border.shp')
print(data)
print(data.keys())
print(data['OBJECTID'])
print(data['BASIN_ID'])
print(data['geometry'])
print("awash crs", data.crs)
#Clip out the composite anomaly using the shapefile
awash_anom = comp_pos_anom.rio.clip(data.geometry.apply(mapping),data.crs)

ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()
ax.add_feature(cartopy.feature.BORDERS,lw=2)
ax.add_feature(cartopy.feature.RIVERS,lw=2)
comp_pos_anom.plot(ax=ax,cmap = plt.cm.BrBG,transform=ccrs.PlateCarree()
#awash_anom.plot(ax=ax,cmap = plt.cm.BrBG,transform=ccrs.PlateCarree()
               #,vmin=-0.8,vmax=0.8
               ,extend='both'
               ,robust=True,cbar_kwargs={'label': "mm/day"})
data.plot(ax=ax, edgecolor='red', facecolor='none',lw=1,zorder=2,linestyle='-')
plt.title(seas+' positive IOD composite rainfall')
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='black', alpha=0.5, linestyle='dotted')
# for extent the order is  [West,East,South,North]
ax.set_extent([31, 48, 3, 15])
gl.top_labels = False
gl.right_labels = False
plt.savefig(path+'plots/C_2.png',bbox_inches='tight',dpi=200)
plt.show()
plt.clf()















