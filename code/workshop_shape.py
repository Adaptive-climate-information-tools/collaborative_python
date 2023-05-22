#Step 1: Start with import statements (this tells python which packages 
#you will actually use in your code and names 
#them, so we will call xarray : xr)
import xarray as xr
from matplotlib import pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature
#from cartopy.io.shapereader import Reader
#from cartopy.feature import ShapelyFeature
from shapely.geometry import mapping
#need to add geopandas and rioxarray to conda environment
import geopandas as gpd
import rioxarray

#%%
#Step 2: Read in our data from a netcdf file of rainfall over Ethiopia 
#(edit the path to the file correctly depending on where you saved it) 
#and select the variable ‘precip’. Now we have a data array called da and
#we can print it to see what we have:
path = '/Users/ellendyer/Library/Mobile Documents/com~apple~CloudDocs/1SHARED_WORK/Work/REACH/Workshop_conda_python/'
ds = xr.open_dataset(path+'workshop_setup/workshop_chirps.nc')
da = ds["precip"]

# print data array
print(da)

#%%
#Step 3 : Xarray has use some of the useful grouping functions like 
#groupby (makes bins), rolling (can do rolling averages), 
#resample (can change time frequencies) and more!

#Let's select months from the JJAS season
da_jjas = da.sel(time=da.time.dt.month.isin([6,7,8,9]))
#Calculate long term mean (of full timeseries for all months)
da_mean = da.mean('time')
#Calculate JJAS mean for full timeseries
da_jjas_mean = da_jjas.mean('time')

#Calculate anomaly of JJAS for each year relative to long term mean
seas_anom = da_jjas_mean - da_mean

#Let's preview this with a plot
seas_anom.plot()
plt.show()
plt.clf()


#%%
#Step 4 : plot the seasonal anomaly map with borders

#But we can make nicer plots with geography (using cartopy)
#We can look at more ways to plot and the details below later on

#---plot with two shapefile approaches (caropy and shapely, and gpd)
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()
ax.add_feature(cartopy.feature.BORDERS)
ax.add_feature(cartopy.feature.RIVERS)
#plot  data
seas_anom.plot(ax=ax,transform=ccrs.PlateCarree())
# for extent the order is  [West,East,South,North]
ax.set_extent([33, 47, 3.5, 14.75])
plt.show()
plt.clf()

#%%

#Plot with a shapefile
#Start by setting up your plot with the standard cartopy method
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()
ax.add_feature(cartopy.feature.BORDERS)
ax.add_feature(cartopy.feature.RIVERS)

#Use geopandas to read in the shapefile (above we imported geopandas as gpd)
data = gpd.read_file(path+'LakeTana_WGS/Lake_Tana_WGS.shp')
print(data.keys())
print("tana crs", data.crs)
data.plot(ax=ax, edgecolor='lightgreen', facecolor='none',lw=2,zorder=2,linestyle='-')
#plot  data
seas_anom.plot(ax=ax,transform=ccrs.PlateCarree())
# for extent the order is  [West,East,South,North]
ax.set_extent([33, 47, 3.5, 14.75])
plt.show()
plt.clf()

#%%
#Tana basin example
#Clip data using shapefiles
#Have to make sure your xarray has a recognised coordinate system 
#that is the same as shapefiles (this one works - in all situations I have found!)
seas_anom_raster = seas_anom.rio.write_crs("epsg:4326")
#You can now clip with the shapefile
data = gpd.read_file(path+'LakeTana_WGS/Lake_Tana_WGS.shp')
seas_anom_clip = seas_anom_raster.rio.clip(data.geometry.apply(mapping),data.crs)
                                 
#Now plot with cartopy
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()
ax.add_feature(cartopy.feature.BORDERS)
ax.add_feature(cartopy.feature.RIVERS)
seas_anom_clip.plot(ax=ax,transform=ccrs.PlateCarree())
# for extent the order is  [West,East,South,North]
ax.set_extent([35, 39, 9, 14.75])
plt.show()
plt.clf()

#%%
#Awash basin example

#This shapefile contains a lot more information so printing a few
#steps helpf you to find what is in the file
data = gpd.read_file(path+'Awash/Awash_basin_border.shp')
print(data)
print(data.keys())
print(data['OBJECTID'])
print(data['BASIN_ID'])
print(data['geometry'])
print("awash crs", data.crs)
#We will clip using the outline of the whole shapefil
seas_anom_clip = seas_anom_raster.rio.clip(data.geometry.apply(mapping),data.crs)
                                 
#Now set up plot using cartopy
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()
ax.add_feature(cartopy.feature.BORDERS)
ax.add_feature(cartopy.feature.RIVERS)

#there isn't a pure outline but you can select different
#regional outlines from the shapefile like this
#data = data.where(data['OBJECTID']==1)
data.plot(ax=ax, edgecolor='darkorange', facecolor='none',lw=2,zorder=2,linestyle='-')

seas_anom_clip.plot(ax=ax,transform=ccrs.PlateCarree())
# for extent the order is  [West,East,South,North]
#ax.set_extent([37.9, 43.5, 7.7, 12.3])
ax.set_extent([33, 47, 3.5, 14.75])
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='black', alpha=0.5, linestyle='dotted')
gl.top_labels = False
gl.right_labels = False


plt.show()
plt.clf()

#%%
#Ethiopia and ICPAC country examples
data = gpd.read_file(path+'afr_g2014_2013_0/afr_g2014_2013_0.shp')
print(data)
print(data.keys())
print(data['geometry'])
print(data['ADM0_NAME'])
print(data['ICPAC'])
#Select the Ethiopian shape to clip with
datanew=data[data['ADM0_NAME']=='Ethiopia']
#datanew=data[data['IGAD']==4]
print(datanew)
print("countries crs", datanew.crs)
#Clip using datanew (where we selected the shape of Ethiopia)
seas_anom_clip = seas_anom_raster.rio.clip(datanew.geometry.apply(mapping),datanew.crs)                             

ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()
ax.add_feature(cartopy.feature.BORDERS)
ax.add_feature(cartopy.feature.RIVERS)
ax.add_feature(cartopy.feature.OCEAN)

#Overlay with all ICPAC country outlines
data_eth = data.where(data['ADM0_NAME']=='Ethiopia')
data_icpac = data.where(data['ICPAC']==14)
data_eth.plot(ax=ax, edgecolor='none', facecolor='beige',lw=2,zorder=2,linestyle='-')
data_icpac.plot(ax=ax, edgecolor='green', facecolor='none',lw=2,zorder=2,linestyle='-')


seas_anom_clip.plot(ax=ax,transform=ccrs.PlateCarree(),add_colorbar=False)
# for extent the order is  [West,East,South,North]
ax.set_extent([20, 90, -5, 17])
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='black', alpha=0.5, linestyle='dotted')
gl.top_labels = False
gl.right_labels = False
ax.set_title('')
plt.show()
plt.clf()



#%%
# write seas_anom to netCDF to save your work 
#for later or share your calculations
#It's always a good idea to use a defined path!
#Based on the last step this will be cut out for Ethiopia
seas_anom_clip.to_netcdf(path+'files/workshop_out_jjas_anom_clip.nc')

