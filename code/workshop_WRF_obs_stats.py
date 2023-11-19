#Step 1: Start with import statements (this tells python which packages 
#you will actually use in your code and names 
#them, so we will call xarray : xr)
import sys
import pandas as pd
import numpy as np
import collections
from re import search
from scipy.stats import pearsonr, spearmanr
from matplotlib import pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature
from sklearn.metrics import mean_squared_error
#%%

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

col_names_to_rename = ["LAT","LON","Year"]
new_col_names = ["lat","lon","Years"]
path = '/Users/ellendyer/Desktop/Nov23_Python/collaborative_python/'

for dataset in ['OBS','WRF']:
    #You can read in spreadsheet station data as a .csv file here:
    st = pd.read_csv(path+'new_files/'+dataset+'_Jun_2021.csv',header=0,index_col=False)
    #strip extra white space around entries
    st = trim_all_columns(st)
    #print head of the dataframe to see column names and values
    #print(st.head())
    #rename columns that usually have spelling mistakes
    st = st.rename(columns=rename_columns(st.columns,col_names_to_rename,new_col_names))
    #print(st.columns)
    #get days in month to use to pivot the spreadsheet using month number from file
    daysinmonth = pd.Timestamp(int(np.nanmean(st.Years.values)), int(np.nanmean(st.Months.values)), 1).days_in_month

    ##Start rearranging the dataframe so we can work with it as an xarray and also
    #so that we can put in datetimes (useful later on)
    #melt moves around the headers so that days are now a column
    newst = pd.melt(st,id_vars=['STATION','lon', 'lat', 'Years','Months'], value_vars=np.arange(1,daysinmonth+1).astype('str'))
    #Can rename columns so that they have generic names
    newst = newst.rename(columns={"variable": "Day", "value": "PRECIP"})
    #print(newst.head())
    #Drop any rows that don't have values (because the month doesn't have 31 days for instance)
    newst = newst.dropna(axis=0)
    #Make a list of dates for all values in the dataframe - could also include the hour
    date = pd.to_datetime(dict(year=newst.Years, month=newst.Months, day=newst.Day), errors='coerce')#,hour=newst.Time))
    #Add this list of datetimes to the dataframe
    newst['time'] = date.values
    #Now we no longer need the columns Years, Month, Day because we made a datetime with them
    newst = newst.drop(columns=['Years','Months','Day','STATION'])
    #Rearrange the columns in the order that makes sense for the work
    newst = newst[['time','lat','lon','PRECIP']]
    #print(newst.head())
    #You can also sort values, here we did it by time, then lat, then lon
    newst.sort_values(by=['time','lat','lon'],inplace=True)
    #print(newst.head())
    #Some values are still 'empty' ie there may have been an entry registered
    #but there was no readable data, so we replace those '' with nans so we
    #can then drop the nans
    newst = newst.replace('', np.nan)
    newst = newst.dropna(axis=0,how='any')
    #Because PRECIP is the variable we want to later plot
    #we want to make sure all the values are marked as floats
    print(dataset,newst.head())
    newst_out = newst.astype({'PRECIP': 'float64'})
    #Now if we want to turn this into an xarray with PRECIP as the variable
    #we can set time, lat, and lon as index values
    newst_toxr_out = newst_out.set_index(['time','lat','lon'])
            
    # ## Create an xarray or output the dataframe to a file
    # #Now the dataframe is setup to be easily converted to an xarray that
    # #we can make a spatial plot of
    xrst_out = newst_toxr_out.to_xarray()['PRECIP']
    #dropping all dimensions with all nan values because xarray has made a mesh
    xrst_out = xrst_out.dropna('lat','all')
    xrst_out = xrst_out.dropna('lon','all')
    #xrst.mean('time', skipna=True).plot()
    #plt.show()
    #plt.clf()
    xrst_out.to_netcdf(path+'files/out_'+dataset+'.nc')
    #We can also output our dataframe in this new
    #organisation to a csv file (other formats available)
    newst_out.to_csv(path+'files/out_'+dataset+'.csv')
    
    if dataset == 'WRF':
        
        wrf_pd = newst_toxr_out
        wrf_xr = xrst_out
        
        
    if dataset == 'OBS':
        obs_pd = newst_toxr_out
        obs_xr = xrst_out

#make sure there are only points in WRF that are also in OBS (there are only a few different points)
#wrf_pd = wrf_pd[wrf_pd.index.isin(obs_pd.index)]
       
#Will use xarray for the statistics because it's a bit easier to loop through spatial dimensions
#First try to drop points where there is no data - ie nans
wrf_stat = wrf_xr.where(obs_xr.notnull(),drop=True)    
obs_stat = obs_xr.where(obs_xr.notnull(),drop=True)    

#Drop points in WRF where there is no observation data because these can't be compared anyway (only a small #)
wrf_stat = wrf_stat.where(obs_stat.notnull(),drop=True) 

#%%
#*****Pearson Correlation Example (correlation at each point with obs)******
#Create numpy arrays to fill with correlations and pvalues at each point and fill with nans for now 
#Use the shape of observations at time 0 to make arrays of the right size
corr_out = np.zeros(obs_stat[0,:,:].shape)*np.nan
pv_out = np.zeros(obs_stat[0,:,:].shape)*np.nan
#loop through all lat lon points and calculate correlation with WRF and OBS timeseries at each point
#enumerate is useful - you can loop through the object (e.g. lat) and use that to select your xarray point
#but you also get an integer counter that starts at 0 to index the numpy array where correlations and pvalues
#are stored

for y,Y in enumerate(obs_stat.lat):
    for x,X in enumerate(obs_stat.lon):
        #using a 'try' statement because some points only have nans (points without stations) and this will skip those
        #and just put a nan value into the correlation and pvalue arrays - this happens in the 'except' below
        try:
            #if there are just a few nan entries (e.g. missing station data) get a list of these indexes 
            #to skip in pearsonr - spearmanr can handle nans but pearsonr can't
            nas = np.logical_or(np.isnan(wrf_stat.sel(lat=Y,lon=X)), np.isnan(obs_stat.sel(lat=Y,lon=X)))
            #calculate correlation and pvalue and put them into the array at the right point
            corr_out[y,x],pv_out[y,x] = pearsonr(obs_stat.sel(lat=Y,lon=X)[~nas],wrf_stat.sel(lat=Y,lon=X)[~nas])
        except:
            corr_out[y,x],pv_out[y,x] = np.nan,np.nan
#To make plotting easy later can put the correlations and pvales from the numpy arrays into
#xarrays using the function copy. Maybe a copy of OBS xarray at the first timestep to copy the 
#metadata like lat and lon and then copy in the data from corr_out and pv_out.
corr_xr = obs_stat[0,:,:].copy(data=corr_out)
pv_xr = obs_stat[0,:,:].copy(data=pv_out)

#****Excercise! What would you do to mask out correlations where the pvalue is greater than 0.1
# where there correlations aren't statisitically significant? How would you plot only statistically
# significant points?
#****

#%%
#*****Root Mean Squared Error Example (at each point with obs)******
#Create numpy arrays to fill with correlations and pvalues at each point and fill with nans for now 
#Use the shape of observations at time 0 to make arrays of the right size
rmse_out = np.zeros(obs_stat[0,:,:].shape)*np.nan
#loop through all lat lon points and calculate correlation with WRF and OBS timeseries at each point
#enumerate is useful - you can loop through the object (e.g. lat) and use that to select your xarray point
#but you also get an integer counter that starts at 0 to index the numpy array where correlations and pvalues
#are stored

for y,Y in enumerate(obs_stat.lat):
    for x,X in enumerate(obs_stat.lon):
        #using a 'try' statement because some points only have nans (points without stations) and this will skip those
        #and just put a nan value into the correlation and pvalue arrays - this happens in the 'except' below
        try:
            #if there are just a few nan entries (e.g. missing station data) get a list of these indexes 
            #to skip in pearsonr - spearmanr can handle nans but pearsonr can't
            nas = np.logical_or(np.isnan(wrf_stat.sel(lat=Y,lon=X)), np.isnan(obs_stat.sel(lat=Y,lon=X)))
            #calculate the RMSE by setting squared=False in the sklearn function mean_squared_error
            rmse_out[y,x] = mean_squared_error(obs_stat.sel(lat=Y,lon=X)[~nas],wrf_stat.sel(lat=Y,lon=X)[~nas],squared=False)
        except:
            rmse_out[y,x] = np.nan
#To make plotting easy later can put the correlations and pvales from the numpy arrays into
#xarrays using the function copy. Maybe a copy of OBS xarray at the first timestep to copy the 
#metadata like lat and lon and then copy in the data from corr_out and pv_out.
rmse_xr = obs_stat[0,:,:].copy(data=rmse_out)

#%%
#*****Mean Bias Example (at each point with obs)******
#Create numpy arrays to fill with correlations and pvalues at each point and fill with nans for now 
#Use the shape of observations at time 0 to make arrays of the right size
bias_out = np.zeros(obs_stat[0,:,:].shape)*np.nan
#loop through all lat lon points and calculate correlation with WRF and OBS timeseries at each point
#enumerate is useful - you can loop through the object (e.g. lat) and use that to select your xarray point
#but you also get an integer counter that starts at 0 to index the numpy array where correlations and pvalues
#are stored

for y,Y in enumerate(obs_stat.lat):
    for x,X in enumerate(obs_stat.lon):
        #using a 'try' statement because some points only have nans (points without stations) and this will skip those
        #and just put a nan value into the correlation and pvalue arrays - this happens in the 'except' below
        try:
            #if there are just a few nan entries (e.g. missing station data) get a list of these indexes 
            #to skip in pearsonr - spearmanr can handle nans but pearsonr can't
            nas = np.logical_or(np.isnan(wrf_stat.sel(lat=Y,lon=X)), np.isnan(obs_stat.sel(lat=Y,lon=X)))
            #calculate the RMSE by setting squared=False in the sklearn function mean_squared_error
            bias_out[y,x] = (wrf_stat.sel(lat=Y,lon=X)[~nas] - obs_stat.sel(lat=Y,lon=X)[~nas]).mean('time')
        except:
            bias_out[y,x] = np.nan
#To make plotting easy later can put the correlations and pvales from the numpy arrays into
#xarrays using the function copy. Maybe a copy of OBS xarray at the first timestep to copy the 
#metadata like lat and lon and then copy in the data from corr_out and pv_out.
bias_xr = obs_stat[0,:,:].copy(data=bias_out)

#%%
#plot the results with xarrays and cartopy
fig = plt.figure(figsize=(5,5))
ax = plt.axes(projection=ccrs.PlateCarree())
#Create a mesh of lat lons using the xarray of station data we made above
#This is necessary because the xarray is a grid instead of a list
lon, lat = np.meshgrid(corr_xr.lon,corr_xr.lat)
#lon, lat = np.meshgrid(rmse_xr.lon,rmse_xr.lat)
#lon, lat = np.meshgrid(bias_xr.lon,bias_xr.lat)

#Use correlation values as the color of scatter points (c) and 
#p-values as size of scatter points (s)
#OR use RMSE values as the collor of scatter points - all the same size
#cmap sets the colormap - using a diverging colormap here so it is
#symmetric around zero because correlations can go from -1 to 1
#set min and max color values to a smaller range of -0.7 to 0.7 so 
#the points are easer to see.
compplot = ax.scatter(x=lon,y=lat,s=80*(1.-pv_xr),c=corr_xr,
                      cmap=plt.cm.PuOr, alpha=0.9,
                      vmin=-.7,vmax=.7,
                      edgecolors='black',transform=ccrs.PlateCarree())

# compplot = ax.scatter(x=lon,y=lat,s=80,c=rmse_xr,
#                       cmap=plt.cm.magma_r, alpha=0.9,
#                       vmin=2,vmax=18,
#                       transform=ccrs.PlateCarree())

# compplot = ax.scatter(x=lon,y=lat,s=80,c=bias_xr,
#                       cmap=plt.cm.RdBu, alpha=0.9,
#                       vmin=-6,vmax=6,
#                       transform=ccrs.PlateCarree())


plt.colorbar(compplot,label='pearson correlation',ax=ax)
#plt.colorbar(compplot,label='RMSE',ax=ax,extend='both')
#plt.colorbar(compplot,label='bias',ax=ax,extend='both')

ax.coastlines()
ax.add_feature(cartopy.feature.BORDERS)
ax.add_feature(cartopy.feature.RIVERS)
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='black', alpha=0.5, linestyle='dotted')
gl.top_labels = False
gl.right_labels = False
#Change these values to plot a different section of the map
ax.set_extent([33.5, 43.5, 3.5, 15])
plt.title('Correlation between OBS and WRF shown with \n color of dots statistical significance shown \n by dot size (bigger=smaller pvalue)')
#plt.title('RMSE between WRF and OBS')
#plt.title('Bias between WRF and OBS')

plt.savefig(path+'plots/WWOS_1.png', bbox_inches='tight',dpi=200)
#plt.savefig(path+'plots/WWOS_2.png', bbox_inches='tight',dpi=200)
#plt.savefig(path+'plots/WWOS_3.png', bbox_inches='tight',dpi=200)

plt.show()
plt.clf()















