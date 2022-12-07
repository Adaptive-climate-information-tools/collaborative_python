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
hsst = hsst.sel(time=slice('1981-01-01','2019-12-31'))
print(hsst)
# convert degree C to Kelvin
hsst = hsst + 273.15

#%%
##Step 3: Select a region and time slice from ts timeseries (you can use date strings
#with xarray). The .sel slice method chooses the nearest point so
#you don't need to know exact bounds

# Ethiopian rainfall
p_eth = da.sel(lat=slice(3,18),lon=slice(33,48))

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

sst_mean = iodi_mean

#%%
##Step 4: calculate correlation of mean SST with Ethiopian rainfall

# Now correlate mean Nino3.4 SSTs with spatial map of Ethiopian rainfall
# Use a numpy function called apply along axis to loop over rainfall
# lat and lon (https://numpy.org/doc/stable/reference/generated/numpy.apply_along_axis.html)
# put 2d array before simple timeseries in apply_along_axis

#Correlation with monthly anomalies
pclim = p_eth.groupby('time.month').mean(dim='time')
p_eth_anom = p_eth.groupby('time.month') - pclim
sstclim = sst_mean.groupby('time.month').mean('time')
sst_mean_anom = sst_mean.groupby('time.month') - sstclim

#Correlation of seasonal anomaly means
p_eth_seas_anom = (p_eth_anom.sel(time=p_eth_anom.time.dt.month.isin([3,4]))).groupby('time.year').mean('time')
sst_mean_seas_anom = (sst_mean_anom.sel(time=sst_mean_anom.time.dt.month.isin([3,4]))).groupby('time.year').mean('time')
correl_map = np.apply_along_axis(spearmanr,0, p_eth_seas_anom,sst_mean_seas_anom)
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
plt.title('Correlation of selected SSTs with Ethiopian rainfall')
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='black', alpha=0.5, linestyle='dotted')
gl.top_labels = False
gl.right_labels = False
plt.savefig(path+'plots/WCN_1.png',bbox_inches='tight',dpi=200)
plt.show()
plt.clf()

#%%
## Step 5: Reverse this and pick a region of Ethiopia (let's do the Awash basin) and correlate with all Indian Ocean SSTs in MAM

#Use monthly anomalies again - can recycle p_eth_anom because it is 2D and cut out with Awash shapefile
#make sure p_eth_anom has specific coordinate system
p_eth_anom.rio.write_crs("epsg:4326", inplace=True)

data = gpd.read_file(path+'Awash/Awash_basin_border.shp')
print(data)
print(data.keys())
print(data['OBJECTID'])
print(data['BASIN_ID'])
print(data['geometry'])
print("awash crs", data.crs)
p_awash_anom = p_eth_anom.rio.clip(data.geometry.apply(mapping),data.crs)
p_awash_anom = p_awash_anom.mean(dim=('lat','lon'))

# have to start again with SSTs because we want them to be 2D now
isst = gsst.sel(latitude=slice(40,-40),longitude=slice(80,180))
isstclim = isst.groupby('time.month').mean('time')
isst_anom = isst.groupby('time.month') - isstclim

#Correlation of seasonal anomaly means for JAS

p_awash_jas = (p_awash_anom.sel(time=p_awash_anom.time.dt.month.isin([7,8,9]))).groupby('time.year').mean('time')
isst_jas = (isst_anom.sel(time=isst_anom.time.dt.month.isin([7,8,9]))).groupby('time.year').mean('time')
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

plt.title('Correlation of selected Ethiopian rainfall with SSTs')
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='black', alpha=0.5, linestyle='dotted')
gl.top_labels = False
gl.right_labels = False
plt.savefig(path+'plots/WCN_2.png',bbox_inches='tight',dpi=200)
plt.show()
plt.clf()






















