#Step 1: Start with import statements (this tells python which packages 
#you will actually use in your code and names 
#them, so we will call xarray : xr)
import xarray as xr
from matplotlib import pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature
import numpy as np
from shapely.geometry import mapping
import pandas as pd
#need to add geopandas and rioxarray to conda environment
import geopandas as gpd
import rioxarray

#Do comparison for the following basin:
sbasin='Awash Awash'
#If you want to examine rainfall in a named basin with streamflow outside
#of the bounds of that region, include the station name here:
sflow_station='Awash below kokadam'
#If you just want to use streamflow stations in the basin comment the line above
#current station names: name = ["Akaki","Awash below kokadam","Awash@Hombole","Berga","Holeta","Mojo@MojoVillage","Mutinicha"]


#%%
#Step 2: Read in our data from a netcdf file of rainfall over Ethiopia 
#(edit the path to the file correctly depending on where you saved it) 
#and select the variable ‘precip’. Now we have a data array called da and
#we can print it to see what we have:
path = '/Users/ellendyer/Library/Mobile Documents/com~apple~CloudDocs/1SHARED_WORK/Work/REACH/Workshop_conda_python/'
da = xr.open_dataset(path+'files/chirps_daily_ethiopia.nc')['precip']
da = da.rio.write_crs("epsg:4326")
emi = xr.open_dataset(path+'files/out_mult_vars_mowie.nc')['PRECIP']
emi = emi.rio.write_crs("epsg:4326")
if sflow_station:
    mowie = pd.read_excel(path+'files/out_flow_mult_mowie.xlsx',header=0,index_col=False)
    mowie = mowie[mowie['station'] == sflow_station]
    mowie = mowie.drop(columns=['station'])
    mowie = mowie.set_index(['time','lat','lon'])
    mowie = mowie.to_xarray()['flow_cumecs']
    mowie = mowie.rio.write_crs("epsg:4326")
else:  
    mowie = xr.open_dataset(path+'files/out_flow_mult_mowie.nc')['flow_cumecs']
    mowie = mowie.rio.write_crs("epsg:4326")

da = da.sel(time=slice('1990-01-01','2010-12-31'),drop=True)
emi = emi.sel(time=slice('1990-01-01','2010-12-31'),drop=True)
mowie = mowie.sel(time=slice('1990-01-01','2010-12-31'),drop=True)


#%%
#Step 3: Readin in shapefiles and select sub-basins
###Contains: Akaki, Awash US Koka, Mojo
pathshape = "mowe_files/Shape_files/"
data = gpd.read_file(path+pathshape+'Sub-basins.shp')
#change the projection so the shapefile has lat lons instead of 
#coordinates like most of the data that will be clipped
data = data.to_crs({'init': 'epsg:4326'}) 
print(data.keys())
print(data['Sub_basin_'])

#Select the Ethiopian basin shape to clip with
sbasin_shape=data[data['Sub_basin_']==sbasin]

#%%
#Step 4: Clip chirps, emi, and mowie data using sub-basin shapefile selection

#Clip using datanew (where we selected the shape of Ethiopia)
da_sbasin = da.rio.clip(sbasin_shape.geometry.apply(mapping),sbasin_shape.crs)                             
emi_sbasin = emi.rio.clip(sbasin_shape.geometry.apply(mapping),sbasin_shape.crs)     
if sflow_station:
    mowie_sbasin = mowie
else:                        
    mowie_sbasin = mowie.rio.clip(sbasin_shape.geometry.apply(mapping),sbasin_shape.crs) 


#%%
#Step 5: compare timeseries of data in the sub-basins

# sub-basin timeseries

# For rainfall you can either do a basin wide mean
#da_sbasin_ts = da_sbasin.mean(dim=('lat','lon'))
# OR a basin wide sum (how much rain accumulates in the whole
# basin at each timestep)
da_sbasin_ts = da_sbasin.sum(dim=('lat','lon'))

mowie_sbasin_ts = mowie_sbasin.mean(dim=('lat','lon'))

#Drop nans in the mowie timeseires
mowie_sbasin_ts = mowie_sbasin_ts.dropna(dim='time')

#Plot the flow and rainfall daily cycles for a few years of data
fig, ax1 = plt.subplots()

color = 'tab:red'
ax1.set_xlabel('time')
ax1.set_ylabel('precip (mm/day)', color=color)
ax1.plot(da_sbasin_ts.time,da_sbasin_ts, color=color)
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:blue'
ax2.set_ylabel('flow (cumecs)', color=color)  # we already handled the x-label with ax1
ax2.scatter(mowie_sbasin_ts.time,mowie_sbasin_ts, color=color,marker='.')
ax2.tick_params(axis='y', labelcolor=color)

ax1.set_xticks(da_sbasin_ts['time'].values[::365])
ax1.set_xticklabels(da_sbasin_ts['time.year'].values[::365], color="k",rotation=45)

ax1.set_title(sbasin+' rainfall/flow comparison')
fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.savefig(path+'plots/MC_1.png', bbox_inches='tight',dpi=200)
plt.show()
plt.clf()

#File write out section
#You can either write out one variable or both to the same file
rain_accu_pd = da_sbasin_ts.to_pandas().to_frame('rain')
#rain_accu_pd.to_excel(path+'files/mc_rain_out.xlsx')
flow_accu_pd = mowie_sbasin_ts.to_pandas().to_frame('flow')
#flow_accu_pd.to_excel(path+'files/mc_flow_out.xlsx')
accu_out_pd = pd.concat([rain_accu_pd,flow_accu_pd],axis=1,join='outer')
accu_out_pd.to_excel(path+'files/mc_combo_out.xlsx')


#%%
#Step 6 annual cumulation comparison plot

#Plot the cummulative sum of flow and rainfall annual cycles for a few years of data

#Can drop values from the chirps timeseries where there is no mowie data
#In example below we keep all timesteps
da_sbasin_ts_m = da_sbasin_ts.where(mowie_sbasin_ts.time==da_sbasin_ts.time,drop=True)

#Use the filtered rainfall to compare with streamflow over the same period
#If you don't want to do this remove the _m in the next line
rain_accu = da_sbasin_ts_m.groupby('time.year').sum('time')
flow_accu = mowie_sbasin_ts.groupby('time.year').sum('time')

fig, ax1 = plt.subplots(figsize=(8,4))
width=0.4
color = 'teal'
ax1.set_xlabel('time')
ax1.set_ylabel('precip (mm/day)', color=color)
ax1.bar(rain_accu.year-width/2,rain_accu,width, color=color,alpha=0.6)
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'purple'
ax2.set_ylabel('flow (cumecs)', color=color)  # we already handled the x-label with ax1
ax2.bar(flow_accu.year+width/2,flow_accu,width, color=color,alpha=0.6)
ax2.tick_params(axis='y', labelcolor=color)

ax1.set_xticks(rain_accu['year'].values)
ax1.set_xticklabels(rain_accu['year'].values, color="k",rotation=45)

ax1.set_title(sbasin+' annual cumulative rainfall/flow comparison')
fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.savefig(path+'plots/MC_2.png', bbox_inches='tight',dpi=200)
plt.show()
plt.clf()

#File write out section
#You can either write out one variable or both to the same file
rain_accu_pd = rain_accu.to_pandas().to_frame('rain')
#rain_accu_pd.to_excel(path+'files/mc_rain_out.xlsx')
flow_accu_pd = flow_accu.to_pandas().to_frame('flow')
#flow_accu_pd.to_excel(path+'files/mc_flow_out.xlsx')
accu_out_pd = pd.concat([rain_accu_pd,flow_accu_pd],axis=1,join='outer')
accu_out_pd.to_excel(path+'files/mc_combo_out.xlsx')


#%%
#Step 7 seasonal cumulation comparison plot

#Plot the cummulative sum of flow and rainfall annual cycles for a few years of data

seas = 'JAS'
mons = [7,8,9]

#create seasonal means
rain_accu_seas = da_sbasin_ts.sel(time=da_sbasin_ts.time.dt.month.isin(mons)).groupby('time.year').sum('time')
flow_accu_seas = mowie_sbasin_ts.sel(time=mowie_sbasin_ts.time.dt.month.isin(mons)).groupby('time.year').sum('time')

fig, ax1 = plt.subplots(figsize=(8,4))
width=0.4
color = 'darkgreen'
ax1.set_xlabel('time')
ax1.set_ylabel('precip (mm/day)', color=color)
ax1.bar(rain_accu_seas.year-width/2,rain_accu_seas,width, color=color,alpha=0.8)
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'brown'
ax2.set_ylabel('flow (cumecs)', color=color)  # we already handled the x-label with ax1
ax2.bar(flow_accu_seas.year+width/2,flow_accu_seas,width, color=color,alpha=0.8)
ax2.tick_params(axis='y', labelcolor=color)

ax1.set_xticks(rain_accu_seas['year'].values)
ax1.set_xticklabels(rain_accu_seas['year'].values, color="k",rotation=45)

ax1.set_title(sbasin+' '+seas+' cumulative rainfall/flow comparison')
fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.savefig(path+'plots/MC_3.png', bbox_inches='tight',dpi=200)
plt.show()
plt.clf()

#File write out section
#You can either write out one variable or both to the same file
rain_accu_pd = rain_accu_seas.to_pandas().to_frame('rain')
#rain_accu_pd.to_excel(path+'files/mc_rain_out.xlsx')
flow_accu_pd = flow_accu_seas.to_pandas().to_frame('flow')
#flow_accu_pd.to_excel(path+'files/mc_flow_out.xlsx')
accu_out_pd = pd.concat([rain_accu_pd,flow_accu_pd],axis=1,join='outer')
accu_out_pd.to_excel(path+'files/mc_combo_out.xlsx')


#%%
#Step 8 multi-seasonal cumulation comparison plot

#Plot the cummulative sum of flow and rainfall annual cycles for a few years of data

#create seasonal mean for rainfall
seas_r = 'MAM'
mons = [3,4,5]
rain_accu_seas = da_sbasin_ts.sel(time=da_sbasin_ts.time.dt.month.isin(mons)).groupby('time.year').sum('time')

#create seasonal mean for streamflow
seas_f = 'JAS'
mons = [7,8,9]
flow_accu_seas = mowie_sbasin_ts.sel(time=mowie_sbasin_ts.time.dt.month.isin(mons)).groupby('time.year').sum('time')

fig, ax1 = plt.subplots(figsize=(8,4))
width=0.4
color = 'darkblue'
ax1.set_xlabel('time')
ax1.set_ylabel(seas_r+' precip (mm/day)', color=color)
ax1.bar(rain_accu_seas.year-width/2,rain_accu_seas,width, color=color,alpha=0.8)
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'coral'
ax2.set_ylabel(seas_f+' flow (cumecs)', color=color)  # we already handled the x-label with ax1
ax2.bar(flow_accu_seas.year+width/2,flow_accu_seas,width, color=color,alpha=0.8)
ax2.tick_params(axis='y', labelcolor=color)

ax1.set_xticks(rain_accu_seas['year'].values)
ax1.set_xticklabels(rain_accu_seas['year'].values, color="k",rotation=45)

ax1.set_title(sbasin+' cumulative rainfall/flow comparison')
fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.savefig(path+'plots/MC_4.png', bbox_inches='tight',dpi=200)
plt.show()
plt.clf()

#File write out section
#You can either write out one variable or both to the same file
rain_accu_pd = rain_accu_seas.to_pandas().to_frame('rain')
#rain_accu_pd.to_excel(path+'files/mc_rain_out.xlsx')
flow_accu_pd = flow_accu_seas.to_pandas().to_frame('flow')
#flow_accu_pd.to_excel(path+'files/mc_flow_out.xlsx')
accu_out_pd = pd.concat([rain_accu_pd,flow_accu_pd],axis=1,join='outer')
accu_out_pd.to_excel(path+'files/mc_combo_out.xlsx')



#%%
#Step 8 compare map plot - show basin for rainfall accumulation and where
#streamflow station is located 

#Select dates in chirps that match with streamflow 
da_sbasin_m = da_sbasin.where(mowie_sbasin.time==da_sbasin.time,drop=True)


fig = plt.figure(figsize=(6,4))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()
ax.add_feature(cartopy.feature.BORDERS)
ax.add_feature(cartopy.feature.RIVERS)
sbasin_shape.plot(ax=ax, edgecolor='black', facecolor='none',lw=1,zorder=2,linestyle='-')
da_sbasin_m.mean('time', skipna=True).plot.pcolormesh(cmap=plt.cm.viridis_r,
                                                   vmin=0,vmax=5,
                                                   cbar_kwargs={'label': "Mean CHIRPS (mm/day)"})

#Create a mesh of lat lons using the xarray of station data we made above
#This is necessary because the xarray is a grid instead of a list
lon, lat = np.meshgrid(mowie_sbasin.lon,mowie_sbasin.lat)
mowie_scat = ax.scatter(x=lon,y=lat,s=60,c=mowie_sbasin.mean('time', skipna=True),
                        vmin=3,vmax=10,
                        cmap=plt.cm.Oranges,edgecolors=None,transform=ccrs.PlateCarree())
plt.colorbar(mowie_scat,label='Flow (cumecs)')
# for extent the order is  [West,East,South,North]
#ax.set_extent([37.5, 40, 7.5, 10])
sbasin_shape.plot(ax=ax, edgecolor='black', facecolor='none',lw=1,zorder=2,linestyle='-')
#seas_anom_clip.plot(ax=ax,transform=ccrs.PlateCarree())
plt.title(sbasin)
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='black', alpha=0.5, linestyle='dotted')
# for extent the order is  [West,East,South,North]
#ax.set_extent([31, 48, 3, 15])
gl.top_labels = False
gl.right_labels = False
fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.savefig(path+'plots/MC_5.png', bbox_inches='tight',dpi=200)
plt.show()
plt.clf()



