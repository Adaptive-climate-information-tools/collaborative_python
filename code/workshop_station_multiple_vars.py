#Step 1: Start with import statements (this tells python which packages 
#you will actually use in your code and names 
#them, so we will call xarray : xr)
import sys
import glob
import pandas as pd
import xarray as xr
import numpy as np
from matplotlib import pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature
import geopandas as gpd
import rioxarray
#from re import search
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
##Step 2: Read in station data from NMA

#Define some important parameters, variables you want to use and columns you want
#to drop or rename for all files here
col_names_to_drop = ['Stations','ID','Monthly total','Elevation']
col_names_to_rename = ["Latitude","Longitude"]
new_col_names = ["lat","lon"]
variables_to_keep = ['PRECIP','RELHUM']
lat1, lat2 = 2,14
lon1, lon2 = 32,50

#Get a list of spreadsheets in the directory - first specify a directory path
dir_p = '/Users/ellendyer/Library/Mobile Documents/com~apple~CloudDocs/1SHARED_WORK/Work/REACH/Workshop_conda_python/station/'
list_f = glob.glob(dir_p+'*.csv')

#You can read in spreadsheet station data as a .csv file here
#make a list to add all the files you read into
dataf_list = []
#loop through the list of files:
for fi in list_f:
    print(fi)
    st = pd.read_csv(fi,header=0,index_col=False)

    #strip extra white space around entries
    st = trim_all_columns(st)
    #print head of the dataframe to see column names and values
    #print(st.head())
    #Select RAINFALL and RELHUM data to use later
    stpr = st[st['Elements'].isin(variables_to_keep)]
    #drop any columns that aren't useful to your task
    #because we are keeping multiple variables don't drop the Elements column here
    stpr = stpr.drop(columns=drop_columns(stpr.columns,col_names_to_drop))
    #rename columns that usually have spelling mistakes
    stpr = stpr.rename(columns=rename_columns(stpr.columns,col_names_to_rename,new_col_names))
    #check the names of columns by printing the keys of the dataframe
    #print(stpr.keys())
    #add your dataframe to the list created before the loop
    dataf_list.append(stpr)
dataf = pd.concat(dataf_list)
print(dataf)
print(dataf.columns)


#%%
##Step 3: Start rearranging the dataframe so we can work with it as an xarray and also
#so that we can put in datetimes (useful later on)
#melt moves around the headers so that days are now a column
newst = pd.melt(dataf,id_vars=['Years', 'lat', 'lon', 'Month', 'Time', 'Elements'], value_vars=np.arange(1,32).astype('str'))
#Can rename columns so that there is a better label for days and for data values
newst = newst.rename(columns={"variable": "Day"})
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
newst = newst[['time','lat','lon','Elements','value']]
#print(newst.head())
#You can also sort values, here we did it by time, then lat, then lon
newst.sort_values(by=['time','lat','lon'],inplace=True)
#Some values are still 'empty' ie there may have been an entry registered
#but there was no readable data, so we replace those '' with nans so we
#can then drop the nans
newst = newst.replace('', np.nan)
#dropna will drop all nan values from the dataframe
newst = newst.dropna(axis=0)
#we want to make sure all the values are marked as floats for math and plotting later
newst = newst.astype({'value': 'float64'})
#select box of region you want data for (sometimes there are random outliers)
newst = newst.loc[(newst.lat>lat1)& (newst.lat<lat2)& (newst.lon>lon1)& (newst.lon<lon2)]

#Now if we want to turn this into an xarray with the chosen Elements (e.g. PRECIP, RELHUM) as variables
#we can set time, lat, and lon as index values
newst_toxr = newst.set_index(['time','lat','lon'])
#Can pivot to have separate columns for different variables now
#for the pandas dataframe we can just set the index to be time
newst = newst.set_index('time')
#remove duplicate values -- this can happen from manually modifying spreadsheets
newst_toxr = newst_toxr[~newst_toxr.index.duplicated()]
#Pivot columns using the variables in Elements column
newst_toxr = newst_toxr.pivot(columns=['Elements'],values='value')
print(newst_toxr.head())



#%%
## Step 4: Create an xarray or output the dataframe to a file
#If you want to keep multiple variables - create an xarray Dataset
xrst = xr.Dataset.from_dataframe(newst_toxr)
#If you only want to keep one variable - create an xarray DataArray and give variable name
#xrst = newst_toxr.to_dataset()['PRECIP']
#dropping all dimensions with all nan values because xarray has made a mesh
xrst = xrst.dropna('lat','all')
xrst = xrst.dropna('lon','all')
#xrst.mean('time', skipna=True).plot()
#plt.show()
#plt.clf()
xrst.to_netcdf('/Users/ellendyer/Library/Mobile Documents/com~apple~CloudDocs/1SHARED_WORK/Work/REACH/Workshop_conda_python/files/out_mult.nc',mode='w')
#We can also output our dataframe in this new
#organisation to a csv file (other formats available)
newst_toxr.to_csv('/Users/ellendyer/Library/Mobile Documents/com~apple~CloudDocs/1SHARED_WORK/Work/REACH/Workshop_conda_python/files/out_mult_vars.csv')


#%%
## Step 5: Select some months from both the CHIRPS and NMA datasets
# Here we will look at PRECIP in March-May
months = [3,4,5]
#select from xarray created above for plotting
xrst = xrst['PRECIP']
#select from dataframe created above for plotting
newst = newst[newst['Elements'] == 'PRECIP']
#Taking long term MAM mean for station data xarray
xrst_mam = xrst.sel(time=np.in1d(xrst['time.month'], months)).mean('time')
#Do the same thing for the dataframe - it is a bit different in format
newst_mam = newst.loc[(newst.index.month.isin([3,4,5]))].groupby(['lat','lon']).mean().reset_index()


#%%
## Step 6: Plot chirps as a background with a station scatter from both the xarray and pandas dataframe

#Create an empty figure - the size is (width,height) and you might have
#to adjust this later
fig = plt.figure(figsize=(8,5))
#Create an axis for your figure and you can use a coordinate projection
#
ax = plt.axes(projection=ccrs.PlateCarree())
#Create a mesh of lat lons using the xarray of station data we made above
#This is necessary because the xarray is a grid instead of a list
lon, lat = np.meshgrid(xrst_mam.lon,xrst_mam.lat)
xrst_p = ax.scatter(x=lon,y=lat,s=60,c=xrst_mam,vmin=0,vmax=4,
                    edgecolors=None,transform=ccrs.PlateCarree())
plt.colorbar(xrst_p,label='MAM Station xarray (mm/day)',ax=ax)


#For a scatter from the dataframe we need to select the index values for lat and lon
newst_p = ax.scatter(x=newst_mam['lon'],
                      y=newst_mam['lat'],
                      s=60,c=newst_mam['value'],vmin=0,vmax=4,
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
ax.set_extent([33, 45, 5, 15.5])
plt.title('Station rainfall scatter in MAM')
plt.savefig('/Users/ellendyer/Library/Mobile Documents/com~apple~CloudDocs/1SHARED_WORK/Work/REACH/Workshop_conda_python/plots/plot_mult1.png', bbox_inches='tight',dpi=200)
plt.show()
plt.clf()


#%%
## Step 6: Add a shapefile to the map

#Read in your shapefile using a library called geopandas - we called it gpd at the start of the code
shape = gpd.read_file('/Users/ellendyer/Library/Mobile Documents/com~apple~CloudDocs/1SHARED_WORK/Work/REACH/Workshop_conda_python/LakeTana_WGS/Lake_Tana_WGS.shp')

fig = plt.figure(figsize=(5,4))
ax = plt.axes(projection=ccrs.PlateCarree())
shape.plot(ax=ax, edgecolor='red', facecolor='none',lw=2,zorder=2,linestyle='--')
#Create a mesh of lat lons using the xarray of station data we made above
#This is necessary because the xarray is a grid instead of a list
lon, lat = np.meshgrid(xrst_mam.lon,xrst_mam.lat)
xrst_p = ax.scatter(x=lon,y=lat,s=60,c=xrst_mam,vmin=0,vmax=4,
                    edgecolors=None,transform=ccrs.PlateCarree())
plt.colorbar(xrst_p,label='MAM Station xarray (mm/day)',ax=ax)
ax.coastlines()
ax.add_feature(cartopy.feature.BORDERS)
ax.add_feature(cartopy.feature.RIVERS)
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='black', alpha=0.5, linestyle='dotted')
gl.top_labels = False
gl.right_labels = False
ax.set_extent([34, 42, 5, 14.5])
plt.title('Station scatter (xarray) \n in MAM with Tana shape')
plt.savefig('/Users/ellendyer/Library/Mobile Documents/com~apple~CloudDocs/1SHARED_WORK/Work/REACH/Workshop_conda_python/plots/plot_mult2.png', bbox_inches='tight',dpi=200)
plt.show()
plt.clf()