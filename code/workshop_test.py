#Step 1: Start with import statements (this tells python which packages 
#you will actually use in your code and names 
#them, so we will call xarray : xr)
import xarray as xr
from matplotlib import pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature



#%%
#Step 2: Read in our data from a netcdf file of rainfall over Ethiopia 
#(edit the path to the file correctly depending on where you saved it) 
#and select the variable ‘precip’. Now we have a data array called da and
#we can print it to see what we have:
path = '/Users/ellendyer/Library/Mobile Documents/com~apple~CloudDocs/1SHARED_WORK/Work/REACH/Workshop_conda_python/'
da = xr.open_dataarray(path+'workshop_setup/workshop_chirps.nc')
# print data array
print(da)
# print attributes in file metadata
print(da.attrs)
# print a list of array dimensions
print(da.dims)
# print a list of array coordinates
print(da.coords)
# print time and lon dim labels
print(da.time)
print(da.lon)


#%%
#Step 3 : Let’s make a quick plot of one month of data. 
#To do this we will select one date using the .sel function 
#and using the coordinate time and its labels. In xarray you can 
#use date strings to select time steps (use the coordinate format 
#that we saw when printing da.time above). If you instead wanted to 
#choose the first time element you would use .isel which lets you select 
#using normal array indices, so the first month would be .isel(time=0).

#use coordinate names
da.sel(time='2018-07-01').plot()
plt.show()
plt.clf()

print(da.sel(lat=34.5,method='nearest'))

#can use time coordinate and datetime to select specific months 
#from the dataset
pseas = da.sel(time=da.time.dt.month.isin([6,7,8,9]))
seasts = pseas.groupby("time.year").mean('time')

sasts_box = seasts.sel(lat=slice(10,15),lon=slice(33,48))
sasts_box_mean = sasts_box.mean('year')
sasts_box_mean.plot()
plt.show()
plt.clf()

#%%
#climatology is called clime
clime = da.sel(time=slice('1991-01-01','2020-12-31'))
clime = clime.groupby('time.month').mean('time')
print(clime)

#year for anomaly = 1998 and is called year_sel
year_sel = da.sel(time=slice('1998-01-01','1998-12-31'))

year_anom = year_sel.groupby('time.month') - clime

year_anom_plot = year_anom.mean('time').plot()
plt.show()
plt.clf()
#%%
#Step 4 : Let’s select and modify data. First we will make a timeseries by taking 
#an average over the dimensions lat and lon. Remember to use the dimension names 
#that we found in printing da.dims above. We can then select a smaller time range 
#for this time series by using .sel again but this time by takin a slice. 
#If you are selecting a single point you either need to know the exact coordinates 
#or you need to ask xarray to select the ‘nearest’ point.

#make a timeseries plot by taking a regional average
ts = da.mean(dim=('lat','lon'))
print(ts)
ts.plot()
plt.show()
plt.clf()

#select a time slice from ts timeseries (you can use date strings
#with xarray)the .sel slice method chooses the nearest point so 
#you don't need to know exact bounds
ts = ts.sel(time=slice('2015-01-01','2017-12-31'))
ts.plot()
plt.show()
plt.clf()

#if you use .sel without slice then you need to say 
#that you want the nearest point
print(da.sel(lat=3.333,lon=40.111,method='nearest'))

#%%
#Step 5 : Xarray has use some of the useful grouping functions like 
#groupby (makes bins), rolling (can do rolling averages), 
#resample (can change time frequencies) and more!

#groupby just makes bins but doesn’t alter your data array  
print(da.groupby("time.season"))

#if you want to alter it you need to act on it and here we do 
#this by taking a time mean with .mean(‘time’)
#calculate long term seasonal means 
#(now 'time' has been replaced by 'season')
print(da.groupby("time.season").mean('time'))
#calculate a timeseries of annual averages 
#(now 'time' has been replaced by 'year')
print(da.groupby("time.year").mean('time'))

# create new data arrays with these time means 
seas_mean = da.groupby("time.season").mean('time')
ann_mean = da.groupby("time.year").mean('time')

#calculate long term seasonal anomalies from the long term annual mean
#you can do simple matrix math to data arrays just like you do with numpy
seas_anom = seas_mean - ann_mean.mean('year')

#%%
#Step 6 : plot the seasonal anomaly maps

#This will plot all four xarray defined seasons (we can talk about 
#selecting different seasons another time!) and plot them with the 
#automatic xarray plotting function
seas_anom.plot(col="season")
plt.show()
plt.clf()

#But we can make nicer plots with geography (using cartopy)
#Here we will just plot for JJA by using .sel and the new coordinate 
#season with the label ‘JJA’
#We can look at more ways to plot and the details below later on
ax = plt.axes(projection=ccrs.PlateCarree())
seas_anom.sel(season='JJA').plot(ax=ax,transform=ccrs.PlateCarree())
ax.coastlines()
ax.add_feature(cartopy.feature.BORDERS)
ax.add_feature(cartopy.feature.RIVERS)
# for extent the order is  [West,East,South,North]
ax.set_extent([33, 47, 3.5, 14.5])
plt.show()
plt.clf()

#%%
#Step 7 : write a data array back to a netcdf file 

# write seas_anom to netCDF to save your work 
#for later or share your calculations
#It's always a good idea to use a defined path!
seas_anom.to_netcdf(path+'files/workshop_out_sa.nc')
