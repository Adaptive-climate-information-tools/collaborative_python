#Step 1: Start with import statements (this tells python which packages 
#you will actually use in your code and names 
#them, so we will call xarray : xr)
import sys
import glob
import os
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
            if jac_sim >= 0.8:
                rename_cols[s1]=cols_newnames[i] 
    print('*****renaming columns: ',rename_cols)
    return rename_cols

def dict_to_dataframe(datain):
    if type(datain)==dict:
        datal=[]
        for k in datain.keys():
            datal.append(datain[k])
        dataout=pd.concat(datal)
    else:
        dataout=datain
    return dataout

def convert_deg_to_decimal(coordinate):
    """
    Convert a lat/lon in degree/minute/sec to decimal

    Parameters
    ----------
    coordinate: string
        format AA:BB:CC D, where D is N,E,S,W

    Returns
    -------
    decimal_coordinate: float
    """
    coordinate = str(int(coordinate))
    print(coordinate)
    if len(coordinate)<6:
        coordinate = '0'+coordinate
    degree,minute,second = coordinate[0:2],coordinate[2:4],coordinate[4:6]
    print(degree,minute,second)
    degree = int(degree)
    minute = int(minute)
    second = int(second)
    decimal_coordinate = (degree+minute/60+second/3600)
    return decimal_coordinate


#%%
##Step 2: Read in station data from NMA

#Define some important parameters, variables you want to use and columns you want
#to drop or rename for all files here
path = '/Users/ellendyer/Library/Mobile Documents/com~apple~CloudDocs/1SHARED_WORK/Work/REACH/Workshop_conda_python/'
col_names_to_drop = ['Stations','ID','Monthly total','Elevation','Gh id','Eg gh id','Eg el abbreviation','Name','Time']
col_names_to_rename = ["Latitude","Longitude","Elements","Years"]
new_col_names = ["lat","lon","Elements","Years"]
variables_to_keep = ['PRECIP']
lat1, lat2 = 2,14
lon1, lon2 = 32,50

#Get a list of spreadsheets in the directory - first specify a directory path
dir_p = 'mowie_files/AWASH climate from AWBA/'
filename = "*.xl*"
list_f = glob.glob(os.path.join(path,dir_p,filename))

#You can read in spreadsheet station data as a .csv file here
#make a list to add all the files you read into
dataf_list = []
#loop through the list of files:
for fi in list_f:
    #print(fi)
    
    #you can read in spreadsheet station data as an excel file here:
    stin = pd.read_excel(fi,header=0,index_col=False,sheet_name=None)
    st = dict_to_dataframe(stin)
    #IF you read in an excel file you need to make sure all the column names are strings
    st.columns = st.columns.map(str)
    #strip extra white space around entries
    st = trim_all_columns(st)
    #print head of the dataframe to see column names and values
    #print(st.head())
    #rename columns that usually have spelling mistakes
    st = st.rename(columns=rename_columns(st.columns,col_names_to_rename,new_col_names))
    #Select RAINFALL and RELHUM data to use later
    stpr = st[st['Elements'].isin(variables_to_keep)]
    #drop any columns that aren't useful to your task
    #because we are keeping multiple variables don't drop the Elements column here
    stpr = stpr.drop(columns=drop_columns(stpr.columns,col_names_to_drop))
    #Modify lat lons so that they are usable
    stpr['lat'] = stpr['lat'].apply(convert_deg_to_decimal)
    stpr['lon'] = stpr['lon'].apply(convert_deg_to_decimal)
    #check the names of columns by printing the keys of the dataframe
    #print(stpr.keys())
    #add your dataframe to the list created before the loop
    print(st.head())
    print(st.columns)
    dataf_list.append(stpr)
dataf = pd.concat(dataf_list)
print(dataf)
print(dataf.columns)



#%%
##Step 3: Start rearranging the dataframe so we can work with it as an xarray and also
#so that we can put in datetimes (useful later on)
#melt moves around the headers so that days are now a column
newst = pd.melt(dataf,id_vars=['Years', 'lat', 'lon', 'Month', 'Elements'], value_vars=np.arange(1,32).astype('str'))
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
#newst = newst.loc[(newst.lat>lat1)& (newst.lat<lat2)& (newst.lon>lon1)& (newst.lon<lon2)]

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
#xrst = xr.Dataset.from_dataframe(newst_toxr)
#If you only want to keep one variable - create an xarray DataArray and give variable name
xrst = newst_toxr.to_xarray()['PRECIP']
#dropping all dimensions with all nan values because xarray has made a mesh
xrst = xrst.dropna('lat','all')
xrst = xrst.dropna('lon','all')
xrst['lat'].attrs['axis'] = 'Y'
xrst['lon'].attrs['axis'] = 'X'
# xrst.mean('time', skipna=True).plot()
# plt.show()
# plt.clf()
xrst.to_netcdf(path+'files/out_mult_vars_mowie.nc',mode='w')
#We can also output our dataframe in this new
#organisation to a csv file (other formats available)
newst_toxr.to_csv(path+'files/out_mult_vars_mowie.csv')

