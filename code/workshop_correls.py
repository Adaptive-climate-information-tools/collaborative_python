#Step 1: Start with import statements (this tells python which packages 
#you will actually use in your code and names 
#them, so we will call xarray : xr)
import sys
import xarray as xr
import pandas as pd
import numpy as np
#Import for Step 3
from matplotlib import pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature
#Import for Step 4
from scipy.stats import spearmanr
#Import for Step 6
import geopandas as gpd
from shapely.geometry import mapping

#%%
##Step 2: Read in our data from a netcdf file of rainfall over Ethiopia
#(edit the path to the file correctly depending on where you saved it) 
#and select the variable ‘precip’. Now we have a data array called da and
#we can print it to see what we have:
da = xr.open_dataarray("../workshop_setup/workshop_chirps.nc")
#convert from mm/month to mm/day
da = da/da.time.dt.daysinmonth
da = da.sel(time=slice('1981-01-01','2019-12-31'))
# print data array
print(da)
hsst = xr.open_dataset("../files/HadISST_sst.nc")['sst']
print(hsst)
# convert degree C to Kelvin
hsst = hsst + 273.15

#%%
##Step 3: Select a region and time slice from ts timeseries (you can use date strings
#with xarray). The .sel slice method chooses the nearest point so
#you don't need to know exact bounds

# Ethiopian rainfall
p_eth = da.sel(lat=slice(3,18),lon=slice(33,48),time=slice('1981-01-01','2019-12-31'))
print(p_eth)
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()
ax.add_feature(cartopy.feature.BORDERS,lw=2)
ax.add_feature(cartopy.feature.RIVERS,lw=2)
p_eth.mean('time').plot()
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='black', alpha=0.5, linestyle='dotted')
gl.top_labels = False
gl.right_labels = False
plt.show()
plt.clf()

# Global SSTs
gsst = hsst.sel(latitude=slice(60,-60),longitude=slice(-180,180),time=slice('1981-01-01','2019-12-31'))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()
ax.add_feature(cartopy.feature.BORDERS,lw=2)
ax.add_feature(cartopy.feature.RIVERS,lw=2)
gsst.mean('time').plot()
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='black', alpha=0.5, linestyle='dotted')
gl.top_labels = False
gl.right_labels = False
plt.show()
plt.clf()

# Western Indian Ocean SSTs
wsst = gsst.sel(latitude=slice(25,0),longitude=slice(45,70))
wsst_mean = wsst.mean(dim=('latitude','longitude'))

# Nino 3.4 SSTs
nsst = gsst.sel(latitude=slice(5,-5),longitude=slice(-170,-120))
nsst_mean = nsst.mean(dim=('latitude','longitude'))
print(nsst)

#%%
##Step 4: calculate correlation of mean SST with Ethiopian rainfall

#Correlation without removing climatology
## What correlation method is best for the variables?
xr.plot.hist(p_eth)
plt.show()
plt.clf()
## Spearman correlation might be best (need to import!)
# Start with mean Ethiopian rainfall and mean Nino3.4 SSTs
p_eth_mean = p_eth.mean(dim=('lat','lon'))
# Use scipy stats spearmanr (https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.spearmanr.html)
correl_value = spearmanr(p_eth_mean, nsst_mean)
print(correl_value)
# Now correlate mean Nino3.4 SSTs with spatial map of Ethiopian rainfall
# Use a numpy function called apply along axis to loop over rainfall
# lat and lon (https://numpy.org/doc/stable/reference/generated/numpy.apply_along_axis.html)
# put 2d array before simple timeseries in apply_along_axis
correl_map = np.apply_along_axis(spearmanr,0, p_eth,nsst_mean)
print(correl_map.shape)
# correl_map isn't an xarray - let's make it one
da_r = xr.DataArray(data=correl_map[0],dims=["lat","lon"],coords={"lat":p_eth.lat,"lon":p_eth.lon})
da_p = xr.DataArray(data=correl_map[1],dims=["lat","lon"],coords={"lat":p_eth.lat,"lon":p_eth.lon})

#Correlation with monthly anomalies
pclim = p_eth.groupby('time.month').mean(dim='time')
p_eth_anom = p_eth.groupby('time.month') - pclim
nsstclim = nsst_mean.groupby('time.month').mean('time')
nsst_mean_anom = nsst_mean.groupby('time.month') - nsstclim
correl_map = np.apply_along_axis(spearmanr,0, p_eth_anom,nsst_mean_anom)
# correl_map isn't an xarray - let's make it one
da_r = xr.DataArray(data=correl_map[0],dims=["lat","lon"],coords={"lat":p_eth.lat,"lon":p_eth.lon})
da_p = xr.DataArray(data=correl_map[1],dims=["lat","lon"],coords={"lat":p_eth.lat,"lon":p_eth.lon})


#Correlation of seasonal anomaly means
p_eth_JJA = (p_eth_anom.sel(time=p_eth_anom.time.dt.month.isin([6,7,8]))).groupby('time.year').mean('time')
nsst_mean_JJA = (nsst_mean_anom.sel(time=nsst_mean_anom.time.dt.month.isin([6,7,8]))).groupby('time.year').mean('time')
correl_map = np.apply_along_axis(spearmanr,0, p_eth_JJA,nsst_mean_JJA)
# correl_map isn't an xarray - let's make it one
da_r = xr.DataArray(data=correl_map[0],dims=["lat","lon"],coords={"lat":p_eth.lat,"lon":p_eth.lon})
da_p = xr.DataArray(data=correl_map[1],dims=["lat","lon"],coords={"lat":p_eth.lat,"lon":p_eth.lon})


#Plot correlation values where pvalues are less than 0.1
da_r_sig = da_r.where(da_p<0.1)
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()
ax.add_feature(cartopy.feature.BORDERS,lw=2)
ax.add_feature(cartopy.feature.RIVERS,lw=2)
da_r_sig.plot(ax=ax,cmap = plt.cm.PuOr,transform=ccrs.PlateCarree()
               ,vmin=-0.8,vmax=0.8
               ,extend='both'
               ,robust=True,cbar_kwargs={'label': ""})
#plt.title('')
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='black', alpha=0.5, linestyle='dotted')
gl.top_labels = False
gl.right_labels = False
#plt.savefig('.png',bbox_inches='tight',dpi=200)
plt.show()
plt.clf()

#%%
## Step 5: Reverse this and pick a region of Ethiopia (let's do the Awash basin) and correlate with all Indian Ocean SSTs in MAM

#Use monthly anomalies again - can recycle p_eth_anom because it is 2D and cut out with Awash shapefile
#make sure p_eth_anom has specific coordinate system
p_eth_anom.rio.write_crs("epsg:4326", inplace=True)

data = gpd.read_file('../Awash/Awash_basin_border.shp')
print(data)
print(data.keys())
print(data['OBJECTID'])
print(data['BASIN_ID'])
print(data['geometry'])
print("awash crs", data.crs)
p_awash_anom = p_eth_anom.rio.clip(data.geometry.apply(mapping),data.crs)
p_awash_anom = p_awash_anom.mean(dim=('lat','lon'))

# have to start again with SSTs because we want them to be 2D now
isst = gsst.sel(latitude=slice(25,-40),longitude=slice(-180,-60))
isstclim = isst.groupby('time.month').mean('time')
isst_anom = isst.groupby('time.month') - isstclim


#Correlation of seasonal anomaly means for JAS

p_awash_jas = (p_awash_anom.sel(time=p_awash_anom.time.dt.month.isin([6,7,8]))).groupby('time.year').mean('time')
isst_jas = (isst_anom.sel(time=isst_anom.time.dt.month.isin([6,7,8]))).groupby('time.year').mean('time')
correl_map = np.apply_along_axis(spearmanr,0, isst_jas,p_awash_jas)
# correl_map isn't an xarray - let's make it one
da_r = xr.DataArray(data=correl_map[0],dims=["lat","lon"],coords={"lat":isst.latitude.values,"lon":isst.longitude.values})
da_p = xr.DataArray(data=correl_map[1],dims=["lat","lon"],coords={"lat":isst.latitude.values,"lon":isst.longitude.values})


#Plot correlation values where pvalues are less than 0.1
da_r_sig = da_r.where(da_p<0.1)

ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()
ax.add_feature(cartopy.feature.BORDERS,lw=2)
da_r.plot(ax=ax,cmap = plt.cm.PuOr,transform=ccrs.PlateCarree()
               ,vmin=-0.8,vmax=0.8
               ,extend='both'
               ,robust=True,cbar_kwargs={'label': ""})
CS = da_p.plot.contour(ax=ax,levels=[0.0,0.05,0.1],transform=ccrs.PlateCarree()
                  ,linestyles=['dotted','--','-']
                  ,colors=['green','black','red'])
ax.clabel(CS, CS.levels, inline=True, fontsize=7)

#plt.title('')
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='black', alpha=0.5, linestyle='dotted')
gl.top_labels = False
gl.right_labels = False
#plt.savefig('.png',bbox_inches='tight',dpi=200)
plt.show()
plt.clf()






















