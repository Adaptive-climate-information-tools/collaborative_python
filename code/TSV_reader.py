"""
TSV_reader.py
Script to extract data from TSV files
"""

import numpy as np
import os
import pandas as pd
from matplotlib import pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature
import sys


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
    data_array = np.array(data_string.strip('\t\n').split('\t')).astype('float')
    return data_array

def process_metadata(metadata, separator_meta='cpt:'):
    """
    Process the line of metadata present in each block of data.
    For now:  only extract the year

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
    metadata = metadata.split(separator_meta)
    date = [s for s in metadata if s.startswith('T=') == True][0]
    date = date.split('=')[1].strip(', \n')
    missing_data = [s for s in metadata if s.startswith('missing=') == True][0]
    missing_data = float(missing_data.split('=')[1].split(',')[0].strip(', \n'))
    units = [s for s in metadata if s.startswith('units=') == True][0]
    units = units.split('=')[1].split(',')[0].strip(', \n')
    time = int(date.split('-')[0])
    return time, missing_data, units

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
        lat.append(line_array[0])
        data_out.append(line_array[1:])
    lat = np.array(lat)
    data_out = np.array(data_out)
    return lat, data_out


def process_block(data_in, index_start, index_stop):
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
    index_start: int
        Indicate where the data block start
    index_stop: int
        Indicate where the data block end (well end+1 - used for slicing)

    Returns
    -------
    data_df:   pandas.DataFrame
        Contains data with longitudes as columns headers, and year and lat as extra columns
    """
    metadata = data_in[index_start]
    year, missing, units = process_metadata(metadata)
    lon = split_line(file_array[index_start + 1])
    data = data_in[index_start + 2:index_stop]
    lat, data = process_data(data)
    data_df = pd.DataFrame(data, columns=lon)
    data_df = data_df.replace(missing,np.nan)
    data_df['lat'] = lat
    data_df['year'] = year
    return data_df, units

#name = 'NASA-GEOSS2S_PRCP_Jul-Sep_iniJun'  # has 4 header rows
#name = 'NextGen_PRCP_Jul-Sep_iniJun'  # has 2 header rows
#name = 'obs_PRCP_Jul-Sep' # has 4 header rows
name = 'NCEP-CFSv2fcst_PRCP_Jul-Sep_iniJun2021' # has 4 header rows

#filepath direct to pycpt directory
#filepath = "/Users/ellendyer/Documents/Work/iri-pycpt/read_tsv/"
#filename = "../PRCP_Jun/input/"+name+".tsv"

#filepath with workshop files:
filepath = "./pycpt_files/"
filename = name+".tsv"
file = os.path.join(filepath, filename)  # File could be entered directly

var = 'PRECIP' # variable name
headerlines = 4  # Defined by the files
block_separator = 'cpt'

header = []
with open(file) as f:
    for i in range(headerlines):
        header.append(f.readline())
    file_array = np.array(f.readlines())

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
    partial_data, units = process_block(file_array, i, i + chan_len)
    full_data_list.append(partial_data)

full_data = pd.concat(full_data_list)
full_data = full_data.melt(id_vars=['year','lat'],var_name='lon',value_name=var)
full_data.set_index(['year', 'lat','lon'], inplace=True)
#print(full_data.head())

#Convert to xarray
full_xr = full_data.to_xarray()[var]
full_xr.attrs['units'] = units
full_xr.attrs['fill_value'] = 'NaN'
print(full_xr)
#sys.exit()

#Write to netcdf
full_xr.to_netcdf('./files/'+name+'.nc')

#Quick plot to check
fig = plt.figure(figsize=(7,6))
ax = plt.axes(projection=ccrs.PlateCarree())
full_xr.sel(year=2021).plot.pcolormesh(cmap=plt.cm.viridis,ax=ax,transform=ccrs.PlateCarree(),
                                alpha=1.0,cbar_kwargs={'label': ""})
ax.coastlines()
ax.add_feature(cartopy.feature.BORDERS)
ax.add_feature(cartopy.feature.RIVERS)
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='black', alpha=0.5, linestyle='dotted')
gl.top_labels = False
gl.right_labels = False
#ax.set_extent([36, 38.5, 9, 13.5])
plt.show()
plt.clf()