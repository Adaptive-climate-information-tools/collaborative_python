"""
MOWIE excel reader
Script to extract data from multiple excel files
"""
import glob
import numpy as np
import os
import pandas as pd
from matplotlib import pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature
#import sys


def split_line(data_string):
    """
    Process data strings from .tsv files

    Parameters
    ----------
    data_string:    string

    Returns
    -------
    data_array:     numpy.array
        Contains the input string broken down and converted to float
    """
    data_array = np.array(data_string.strip('\t\n').split(','))
    return data_array
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
    direction = coordinate[-1]
    coordinate = coordinate[:-1]
    sign = 1
    if direction == 'W' or direction == 'S':
        sign = -1

    degree,minute,second = coordinate.split(':')
    degree = int(degree)
    minute = int(minute)
    second = int(second)
    decimal_coordinate = sign*(degree+minute/60+second/3600)
    return decimal_coordinate
def process_metadata(metadata,year_marker = 'Station Number', lat_lon_marker ='Latitude'):
    """
    Process the metadata
    For now, only year

    Parameters
    ----------
    metadata:   string
    separator_meta: string, optional
        string used to separate the different field. Available in case variable between files
    Returns
    -------
    year:  integer
    missing_data:   float
        value for missing data
    units:  string
    """
    
    index_year = np.where(np.char.startswith(metadata, year_marker))[0][0]
    year_list = metadata[index_year]
    yearindex = np.where(np.char.startswith(year_list,'Year:'))[0][0]
    yearstring = year_list[yearindex]
    year = int(yearstring.split('Year:')[1].split(',')[0])
    
    index_lat_lon = np.where(np.char.startswith(metadata, lat_lon_marker))[0][0]
    lat_lon_list = metadata[index_lat_lon]
    latindex = np.where(np.char.startswith(lat_lon_list,'Latitude'))[0][0]
    latstring = lat_lon_list[latindex]
    lat = latstring.split('Latitude :')[1].split(',')[0].strip()
    lat = convert_deg_to_decimal(lat)
  
    lonindex = np.where(np.char.startswith(lat_lon_list,'Longitude'))[0][0]
    lonstring = lat_lon_list[lonindex]
    lon = lonstring.split('Longitude :')[1].split(',')[0].strip()
    lon = convert_deg_to_decimal(lon)

    return year, lat, lon

def process_data(data_in):
    """
     Process block of data.
     Line by line, converts the data string to a numpy.array of floats
     Assign first value to a latitude array, and the rest to a data array

     Parameters
     ----------
     data_in:   numpy.array
        contains a series of lat + data stringa

     Returns
     -------
     lat:   numpy.array
     data_out:     numpy.ndarray (2d)
        Contains data (float), row representing diff latitudes, columns, longitudes
     """
    lat = []
    data_out = []
    for line in data_in:
        line_array = split_line(line)
        lat.append(int(line_array[0]))
        data_temp = line_array[2:]
        data_float = np.ones(len(data_temp))*np.nan
        ind_valid = np.where((~np.char.isspace(data_temp))&(data_temp!='')&(data_temp!='-')) #Question do we differientiate betweem missing data and out of date slot
        data_float[ind_valid] = data_temp[ind_valid].astype('float')
        # ind_no_day = np.where(np.char.isspace(data_temp))
        # data_float[ind_no_day] = -999

        data_out.append(data_float)
    lat = np.array(lat)
    data_out = np.array(data_out)
    return lat, data_out


def process_block(data_in, block_start, block_stop, filename):
    """
    For each block of data, transform the data into a pandas dataframe
    by extracting:
        - the year from the metadata
        - longitudes from first row of data
        - latitudes from first column of data
        - data from the rest

    Parameters
    ----------
    data_in:   numpy.array
        contains strings for a block of data
    block_start: int
        Indicate where the data block start
    block_stop: int
        Indicate where the data block end (well end+1 - used for slicing)
    filename:  string
        filepath of the whole data - needed to read in the data
    Returns
    -------
    data_df:   pandas.DataFrame
        Contains data with longitudes as columns headers, and year and lat as extra columns
    """
    # The metadata is all around the block - will need to decide what we do
    metadata = data_in[block_start:block_stop]
    year, lat, lon =  process_metadata(metadata)

    local_start = np.where(np.char.startswith(data_in[block_start:block_stop], '1'))[0][0]
    local_stop = np.where(np.char.startswith(data_in[block_start:block_stop], '31'))[0][0]

    global_start = block_start + local_start + 4 - 2
    num_rows = local_stop -  local_start + 2
    #LATEST IDEA!!!
    df= pd.read_excel(filename,skiprows =global_start,nrows=num_rows)
    df = df.drop(columns ='Unnamed: 1')
    df = df.drop(index=0)
    # new_columns ={'Unnamed: 0':'day','Jan':1, 'Feb':2, 'Mar':3, 'Apr':4,'May':5,'Jun':6,'Jul':7,'Aug8,9,10,11,12]
    df.columns = ['Day',1,2,3,4,5,6,7,8,9,10,11,12]
    df['Years'] = year
    df['lat'] = lat
    df['lon'] = lon
    data_df = pd.melt(df, id_vars=['Years', 'lat', 'lon', 'Day'],
                    value_vars=np.arange(1, 13))
    data_df = data_df.rename(columns={"variable": "Month","value":var})
    #
    loc_ = np.where(data_df[var].astype('str').str.isspace())[0]
    data_df = data_df.drop(loc_)
    loc_ = np.where(data_df[var].astype('str').str.contains('-'))[0]
    # data_df.value[loc_] = np.nan# np.ones(len(loc_))*np.nan
    date = pd.to_datetime(dict(year=data_df.Years, month=data_df.Month, day=data_df.Day), errors='coerce')
    data_df['time'] = date.values
    data_df = data_df.drop(columns=['Years', 'Month', 'Day'])
    data_df = data_df[['time', 'lat', 'lon', var]]
    return data_df

#These files are missing lat lon lines/values: "Awash u.s. kokadam","Awash@Ginchi","Teji"
name = ["Akaki","Awash below kokadam","Awash@Hombole","Berga","Holeta","Mojo@MojoVillage","Mutinicha"]
path = '/Users/ellendyer/Library/Mobile Documents/com~apple~CloudDocs/1SHARED_WORK/Work/REACH/Workshop_conda_python/'
filepath = "mowie_files/"
###END OF SECTION OF CODE TO EDIT BEFORE YOU RUN###


var = 'flow_cumecs' # variable name
headerlines = 4  # Defined by the files
block_separator = 'Station Number '

#You can read in spreadsheet station data as a .csv file here
#make a list to add all the files you read into
dataf_list = []
#loop through the list of files:
for N in name:
    filename = N+".xl*"
    file = glob.glob(os.path.join(path,filepath, filename))[0]  # File could be entered directly
    #print(file)

    file_excel = pd.read_excel(file,skiprows=headerlines)
    file_array = file_excel.to_numpy()
    file_array = file_array.astype(str)
    
    # Find the beginning of each block by look for the block separator
    indices = np.where(np.char.startswith(file_array, block_separator))[0]
    # Define the length of each block - it assumes all blocks are the same length
    if len(indices)>1:
        chan_len = np.diff(indices)[0]
    else:
        chan_len = len(file_array)
    
    num_block = len(indices)
    full_data_list = []
    
    for n in range(num_block):
        """
        Process a block at a time
        Put data dataframe in a list
        """
        i = indices[n]
        partial_data = process_block(file_array, i, i + chan_len,file)
        partial_data['station'] = N
        #print('**********', partial_data)
        full_data_list.append(partial_data)
    
    full_data = pd.concat(full_data_list)
    full_data = full_data.replace({'-   ':np.nan}) #Hack solution for now
    full_data = full_data.dropna(subset='time')
    print(N,full_data.head(1))
    print('**********************************')

    dataf_list.append(full_data)
dataf = pd.concat(dataf_list)


#Preparing data to be exported as an xarray or a pandas data frame

#Now if we want to turn this into an xarray with PRECIP as the variable
#we can set time, lat, and lon as index values
data_toxr = dataf.set_index(['time','lat','lon'])
#for the pandas dataframe we can just set the index to be time
data_topd = dataf.set_index('time')

#Now the dataframe is setup to be easily converted to an xarray that
#we can make a spatial plot of
xrst = data_toxr.to_xarray()[var]
#dropping all dimensions with all nan values because xarray has made a mesh
xrst = xrst.dropna('lat','all')
xrst = xrst.dropna('lon','all')
xrst['lat'].attrs['axis'] = 'Y'
xrst['lon'].attrs['axis'] = 'X'
xrst.mean('time', skipna=True).plot()
plt.show()
plt.clf()
xrst.to_netcdf(path+'files/out_flow_mult_mowie.nc',mode='w')
#We can also output our dataframe in this new
#organisation to an excel file (other formats available)
data_topd.to_excel(path+'files/out_flow_mult_mowie.xlsx')


