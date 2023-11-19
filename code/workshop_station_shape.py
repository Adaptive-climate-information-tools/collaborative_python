#Step 1: Start with import statements (this tells python which packages 
#you will actually use in your code and names 
#them, so we will call xarray : xr)
#import sys
import xarray as xr
import pandas as pd
import numpy as np
#Import for Step 7
from matplotlib import pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature
#Import for Step 9
#**To run this code you will need to conda install geopandas and rioxarray
import geopandas as gpd
from shapely.geometry import mapping
import rioxarray
import collections

#Here is a function definition - someone wrote this and shared it
#to help other strip empty spaces from excel file fields
#and we will 'call it' later in the code. Usually function definitions
#are put at the start of the code.
def trim_all_columns(df):
    """
    Trim whitespace from ends of each value across all series in dataframe
    """
    trim_strings = lambda x: x.strip() if isinstance(x, str) else x
    return df.applymap(trim_strings)

def drop_columns(cols_indf,cols_todrop):
    """Computes the Jaccard similarity score between s1 and s2.
    https://mindee.com/blog/partial-string-matching/
    https://en.wikipedia.org/wiki/Jaccard_index
    """    
    drop_cols = []
    for s1 in cols_indf:
        for s2 in cols_todrop:
            jac_sim = len(set(s1.lower()) & set(s2.lower())) / len(set(s1.lower()) | set(s2.lower()))
            #print(s1,s2,jac_sim)
            if jac_sim > 0.8:
                drop_cols.append(s1) 
    print('*****dropping columns: ',drop_cols)
    return drop_cols

def rename_columns(cols_indf,cols_torename,cols_newnames):
    """Same Jaccard method as drop_columns"""
    rename_cols = collections.OrderedDict()
    for s1 in cols_indf:
        for i,s2 in enumerate(cols_torename):
            jac_sim = len(set(s1.lower()) & set(s2.lower())) / len(set(s1.lower()) | set(s2.lower()))
            #print(s1,s2,jac_sim)
            if jac_sim > 0.9:
                rename_cols[s1]=cols_newnames[i] 
    print('*****renaming columns: ',rename_cols)
    return rename_cols

#%%
##Step 2: Read in our data from a netcdf file of rainfall over Ethiopia
#(edit the path to the file correctly depending on where you saved it) 
#and select the variable ‘precip’. Now we have a data array called da and
#we can print it to see what we have:
path = '/Users/ellendyer/Desktop/Nov23_Python/collaborative_python/'
da = xr.open_dataarray(path+'workshop_setup/workshop_chirps.nc')
#convert from mm/month to mm/day
da = da/da.time.dt.daysinmonth
da = da.sel(time=slice('1996-01-01','2010-12-31'))
# print data array
#print(da)

#select a region and time slice from ts timeseries (you can use date strings
#with xarray). The .sel slice method chooses the nearest point so
#you don't need to know exact bounds
chirps_reg = da.sel(lat=slice(7,12),lon=slice(36,40))
#print(chirps_reg)

#%%
##Step 3: Read in station data from NMA

col_names_to_drop = ['Stations','ID','Monthly total','Elevation','Elements']
col_names_to_rename = ["Latitude","Longitude"]
new_col_names = ["lat","lon"]
variables_to_keep = 'PRECIP'
#You can read in spreadsheet station data as a .csv file here:
#st = pd.read_csv(path+'station/NMA_Tana_basin.csv',header=0,index_col=False)

#OR you can read in spreadsheet station data as an excel file here:
path = '/Users/ellendyer/Desktop/Nov23_Python/collaborative_python/'
st = pd.read_excel(path+"station/NMA_Tana_basin.xlsx",header=0,index_col=False,sheet_name='data')
#IF you read in an excel file you need to make sure all the column names are strings
st.columns = st.columns.map(str)

#strip extra white space around entries
st = trim_all_columns(st)
#print head of the dataframe to see column names and values
print(st.head())
#Select just rainfall data (here called 'PRECIP')
stpr = st[st['Elements'] == variables_to_keep]
#drop any columns that aren't useful to your task
stpr = stpr.drop(columns=drop_columns(stpr.columns,col_names_to_drop))
#rename columns that usually have spelling mistakes
stpr = stpr.rename(columns=rename_columns(stpr.columns,col_names_to_rename,new_col_names))

#check the names of columns by printing the keys of the dataframe
print(stpr.keys())

#%%
##Step 4: Start rearranging the dataframe so we can work with it as an xarray and also
#so that we can put in datetimes (useful later on)
#melt moves around the headers so that days are now a column
newst = pd.melt(stpr,id_vars=['Years', 'lat', 'lon', 'Month', 'Time'], value_vars=np.arange(1,32).astype('str'))
#Can rename columns so that they have generic names
newst = newst.rename(columns={"variable": "Day", "value": variables_to_keep})
print(newst.head())
#Drop any rows that don't have values (because the month doesn't have 31 days for instance)
newst = newst.dropna(axis=0)
#Make a list of dates for all values in the dataframe - could also include the hour
date = pd.to_datetime(dict(year=newst.Years, month=newst.Month, day=newst.Day), errors='coerce')#,hour=newst.Time))
#Add this list of datetimes to the dataframe
newst['time'] = date.values
#Now we no longer need the columns Years, Month, Day because we made a datetime with them
newst = newst.drop(columns=['Years','Month','Day'])
#Rearrange the columns in the order that makes sense for the work
newst = newst[['time','lat','lon',variables_to_keep]]
#print(newst.head())
#You can also sort values, here we did it by time, then lat, then lon
newst.sort_values(by=['time','lat','lon'],inplace=True)
#Some values are still 'empty' ie there may have been an entry registered
#but there was no readable data, so we replace those '' with nans so we
#can then drop the nans
newst = newst.replace('', np.nan)
newst = newst.dropna(axis=0)
#Because PRECIP is the variable we want to later plot
#we want to make sure all the values are marked as floats
newst = newst.astype({variables_to_keep: 'float64'})
#Now if we want to turn this into an xarray with PRECIP as the variable
#we can set time, lat, and lon as index values
newst_toxr = newst.set_index(['time','lat','lon'])
#Select same date range as chirps
newst_toxr = newst_toxr.loc['1996-01-01':'2010-12-31']
print(newst_toxr.head())

#%%
## Step 5: Create an xarray or output the dataframe to a file
#Now the dataframe is setup to be easily converted to an xarray that
#we can make a spatial plot of
xrst = newst_toxr.to_xarray()[variables_to_keep]
#dropping all dimensions with all nan values because xarray has made a mesh
xrst = xrst.dropna('lat','all')
xrst = xrst.dropna('lon','all')
#xrst.mean('time', skipna=True).plot()
#plt.show()
#plt.clf()
xrst.to_netcdf(path+'files/out.nc')
#We can also output our dataframe in this new
#organisation to a csv file (other formats available)
newst.to_csv(path+'files/out.csv')


#%%
## Step 6 : Select some months from both the CHIRPS and NMA datasets
# Here we will look at March-May
months = [3,4,5]
#Taking long term MAM mean for chirps
da_mam = da.sel(time=np.in1d(da['time.month'], months)).mean('time')

#Taking long term MAM mean for station data xarray
xrst_mam = xrst.sel(time=np.in1d(xrst['time.month'], months)).mean('time')
newst_mam = newst_toxr.loc[(newst_toxr.index.get_level_values('time').month.isin([3,4,5]))].groupby(level=('lat','lon')).mean()
#print(newst_mam.index.get_level_values('lat').values)
#print(newst_mam.head(15))

#%%
## Step 7: Plot chirps as a background with a station scatter from both the xarray and pandas dataframe

#Create an empty figure - the size is (width,height) and you might have
#to adjust this later
fig = plt.figure(figsize=(8,5))
#Create an axis for your figure and you can use a coordinate projection
#
ax = plt.axes(projection=ccrs.PlateCarree())
#Do a filled contour with chirps (same as first example)
#alpha of less than 1 makes the contours a bit transparent so its easier to see the scatter
da_mam.plot.contourf(ax=ax,transform=ccrs.PlateCarree(),alpha=0.8,vmin=0,vmax=4,
                     cbar_kwargs={'label': "MAM CHIRPS contours (mm/day)"})
#Create a mesh of lat lons using the xarray of station data we made above
#This is necessary because the xarray is a grid instead of a list
lon, lat = np.meshgrid(xrst_mam.lon,xrst_mam.lat)
xrst_p = ax.scatter(x=lon,y=lat,s=60,c=xrst_mam,vmin=0,vmax=4,
                    edgecolors=None,transform=ccrs.PlateCarree())
plt.colorbar(xrst_p,label='MAM Station xarray (mm/day)',ax=ax)
#For a scatter from the dataframe we need to select the index values for lat and lon
#Because there are multiple indices (time, lat, lon) we need to use get_level_values
#because these indices are called values
newst_p = ax.scatter(x=newst_mam.index.get_level_values('lon').values,
                     y=newst_mam.index.get_level_values('lat').values,
                     s=60,c=newst_mam[variables_to_keep].values,vmin=0,vmax=4,
                     edgecolors='black',transform=ccrs.PlateCarree())
plt.colorbar(newst_p,label='MAM Station dataframe (mm/day)',ax=ax)
#*** You can't see both scatters because they are on top of each other - comment each one
#    out at a time and plot to see if you can see a difference in them
#Now we add all the cartopy information we did last time and add in some gridlines so it
#is easier to read the plot
ax.coastlines()
ax.add_feature(cartopy.feature.BORDERS)
ax.add_feature(cartopy.feature.RIVERS)
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='black', alpha=0.5, linestyle='dotted')
gl.top_labels = False
gl.right_labels = False
# for extent the order is  [West,East,South,North]
ax.set_extent([36, 38.5, 9, 13.5])
plt.title('Station scatter (pandas and xarray) \n over CHIRPS rainfall contours \n in MAM')
plt.savefig(path+'plots/WSS_1.png', bbox_inches='tight',dpi=200)
plt.show()
plt.clf()

#%%
## Step 8: Interpolate Chirps data to the station locations so we can then calculate the difference
#between values at station locations
#There are options for interpolation methods but I have chosen to linearly interpolate here
da_st = da.interp(lon=xrst.lon.values,lat=xrst.lat.values,method='linear')
diff_st = da_st - xrst
diff_st = diff_st.mean('time')
#print(diff_st)

fig = plt.figure(figsize=(7,6))
ax = plt.axes(projection=ccrs.PlateCarree())
#Do a filled contour with chirps (same as first example)
#alpha of less than 1 makes the contours a bit transparent so its easier to see the scatter
da.mean('time').plot.pcolormesh(cmap=plt.cm.viridis,ax=ax,transform=ccrs.PlateCarree(),
                                alpha=1.0,cbar_kwargs={'label': "Mean CHIRPS mesh (mm/day)"})
#Create a mesh of lat lons using the xarray of station data we made above
#This is necessary because the xarray is a grid instead of a list
lon, lat = np.meshgrid(diff_st.lon,diff_st.lat)
diffp = ax.scatter(x=lon,y=lat,s=60,c=diff_st,vmin=-1.5,vmax=1.5,
                   cmap=plt.cm.PRGn,edgecolors=None,transform=ccrs.PlateCarree())
plt.colorbar(diffp,label='Mean CHIRPS-station (mm/day)',ax=ax)
ax.coastlines()
ax.add_feature(cartopy.feature.BORDERS)
ax.add_feature(cartopy.feature.RIVERS)
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='black', alpha=0.5, linestyle='dotted')
gl.top_labels = False
gl.right_labels = False
ax.set_extent([36, 38.5, 9, 13.5])
plt.title('Difference between interpoloated \n CHIRPS and station '
          'data \n over CHIRPS rainfall mesh \n for long term mean')
plt.savefig(path+'plots/WSS_2.png', bbox_inches='tight',dpi=200)
plt.show()
plt.clf()


#%%
## Step 9: Try to cut out this data for Tana basin where the NMA station data is located

#Read in your shapefile using a library called geopandas - we called it gpd at the start of the code
shape = gpd.read_file(path+'LakeTana_WGS/Lake_Tana_WGS.shp')
#Now you have to add some coordinate information to your CHIRPS datarray (in this case epsg:4326)
da_rio = da.rio.write_crs("epsg:4326", inplace=True)
#Now that the shapefile and dataframe have coordinates you can use rio.clip to cut the dataarray out in
#using the shapefile
da_clip = da.rio.clip(shape.geometry.apply(mapping), shape.crs)

fig = plt.figure(figsize=(7,6))
ax = plt.axes(projection=ccrs.PlateCarree())
#Do a filled contour with chirps (same as first example)
#alpha of less than 1 makes the contours a bit transparent so its easier to see the scatter
da_clip.mean('time').plot.pcolormesh(cmap=plt.cm.viridis,ax=ax,transform=ccrs.PlateCarree(),
                                alpha=1.0,cbar_kwargs={'label': "Mean CHIRPS mesh (mm/day)"})
shape.plot(ax=ax, edgecolor='black', facecolor='none',lw=2,zorder=2,linestyle='--')
#Create a mesh of lat lons using the xarray of station data we made above
#This is necessary because the xarray is a grid instead of a list
lon, lat = np.meshgrid(diff_st.lon,diff_st.lat)
diffp = ax.scatter(x=lon,y=lat,s=60,c=diff_st,vmin=-1.5,vmax=1.5,
                   cmap=plt.cm.PRGn,edgecolors=None,transform=ccrs.PlateCarree())
plt.colorbar(diffp,label='Mean CHIRPS-station (mm/day)')
ax.coastlines()
ax.add_feature(cartopy.feature.BORDERS)
ax.add_feature(cartopy.feature.RIVERS)
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='black', alpha=0.5, linestyle='dotted')
gl.top_labels = False
gl.right_labels = False
ax.set_extent([36, 38.5, 9, 13.5])
plt.title('Difference between interpoloated \n CHIRPS and station '
          'data \n over CHIRPS rainfall mesh \n for long term mean')
plt.savefig(path+'plots/WSS_3.png', bbox_inches='tight',dpi=200)
plt.show()
plt.clf()