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
import numpy as np

#%%
#Step 2: Read in our data from a netcdf file of rainfall over Ethiopia 
#(edit the path to the file correctly depending on where you saved it) 
#and select the variable ‘precip’. Now we have a data array called da and
#we can print it to see what we have:
path = '/Users/ellendyer/Library/Mobile Documents/com~apple~CloudDocs/1SHARED_WORK/Work/REACH/Workshop_conda_python/'
ds = xr.open_dataset(path+'files/chirps_monthly_ethiopia.nc')
da = ds["precip"].mean('time')
da = da.rio.write_crs("epsg:4326")
# print data array
print(da)

pathshape = "mowie_files/Shape_files/"
data = gpd.read_file(path+pathshape+'Awash_River_Basin.shp')
print(data.keys())
print(data['Name'].values)
datanew=data[data['Name']=='Awash']
data_crs = data.to_crs({'init': 'epsg:4326'}) 
fig = plt.figure(figsize=(14,7))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()
ax.add_feature(cartopy.feature.BORDERS)
ax.add_feature(cartopy.feature.RIVERS)
da.plot.pcolormesh(cmap=plt.cm.Blues,ax=ax,transform=ccrs.PlateCarree(),vmin=0,vmax=3,
                                alpha=0.5,cbar_kwargs={'label': "Mean CHIRPS (mm/day)"})
data_crs.plot(ax=ax, edgecolor='black', facecolor='none',lw=1,zorder=2,linestyle='-')
emivar = xr.open_dataset(path+'files/out_mult_vars_mowie.nc')
mowievar = xr.open_dataset(path+'files/out_flow_mult_mowie.nc')
emivar['PRECIP'].mean('time', skipna=True).plot(cmap=plt.cm.viridis_r,vmin=0,vmax=3)
#mowievar['flow_cumecs'].mean('time', skipna=True).plot(cmap=plt.cm.Reds,vmin=0,vmax=35)

#Create a mesh of lat lons using the xarray of station data we made above
#This is necessary because the xarray is a grid instead of a list
lon, lat = np.meshgrid(mowievar['flow_cumecs'].lon,mowievar['flow_cumecs'].lat)
mowie_scat = ax.scatter(x=lon,y=lat,s=60,c=mowievar['flow_cumecs'].mean('time', skipna=True),
                        vmin=0,vmax=35,
                        cmap=plt.cm.Reds,edgecolors=None,transform=ccrs.PlateCarree())
ax.annotate('Awash below Koka',xy=(39.2,8.47),xytext=(39.25,8.47),fontsize='x-small',zorder=3)
ax.annotate('Awash@Hombole',xy=(38.8,8.38),xytext=(38.85,8.38),fontsize='x-small',zorder=3)
plt.colorbar(mowie_scat,label='Flow (cumecs)')




# for extent the order is  [West,East,South,North]
#ax.set_extent([36.5, 43.5, 7.5, 12.5])
ax.set_extent([37.5, 41, 7.5, 10])

#Step 3: Read in shapefile and clip data with the basin or region of choice

###Contains: Akaki, Awash US Koka, Mojo
#pathshape = "mowie_files/shape file awash awash sub basin/"
#data = gpd.read_file(path+pathshape+'basin_awash.shp')
###Contains: Akaki, Awash US Koka, Mojo
pathshape = "mowie_files/Shape_files/"
data = gpd.read_file(path+pathshape+'Sub-basins.shp')
#change the projection so the shapefile has lat lons instead of 
#coordinates like most of the data that will be clipped
data = data.to_crs({'init': 'epsg:4326'}) 
print(data.keys())
print(data['Sub_basin_'])


#Select the Ethiopian basin shape to clip with
Akaki_shape=data[data['Sub_basin_']=='Akaki']
Mojo_shape=data[data['Sub_basin_']=='Mojo']
USKoka_shape=data[data['Sub_basin_']=='Awash US Koka']
AwashAwash_shape=data[data['Sub_basin_']=='Awash Awash']
#Clip using datanew (where we selected the shape of Ethiopia)
#seas_anom_clip = da.rio.clip(datanew_crs.geometry.apply(mapping),datanew_crs.crs)                             

Akaki_shape.plot(ax=ax, edgecolor='red', facecolor='none',lw=1,zorder=2,linestyle='-')
Mojo_shape.plot(ax=ax, edgecolor='red', facecolor='none',lw=1,zorder=2,linestyle='-')
USKoka_shape.plot(ax=ax, edgecolor='red', facecolor='none',lw=1,zorder=2,linestyle='-')
AwashAwash_shape.plot(ax=ax, edgecolor='red', facecolor='none',lw=1,zorder=2,linestyle='-')
#seas_anom_clip.plot(ax=ax,transform=ccrs.PlateCarree())
#plt.savefig(path+'plots/MS_1.png', bbox_inches='tight',dpi=200)

gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='black', alpha=0.5, linestyle='dotted')
# for extent the order is  [West,East,South,North]
#ax.set_extent([31, 48, 3, 15])
gl.top_labels = False
gl.right_labels = False


plt.savefig(path+'plots/MS_2.png', bbox_inches='tight',dpi=200)
plt.show()
plt.clf()

