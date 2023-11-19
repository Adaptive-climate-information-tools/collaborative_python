#Step 1: Start with import statements (this tells python which packages 
#you will actually use in your code and names 
#them, so we will call xarray : xr)
import sys
import xarray as xr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#%%
##Step 2: Read in our data from a netcdf file of rainfall over Ethiopia
#(edit the path to the file correctly depending on where you saved it) 
#and select the variable ‘precip’. Now we have a data array called da and
#we can print it to see what we have:
path = '/Users/ellendyer/Desktop/Nov23_Python/collaborative_python/'
da = xr.open_dataarray(path+'workshop_setup/workshop_chirps.nc')
#convert from mm/month to mm/day
da = da/da.time.dt.daysinmonth
da = da.sel(time=slice('1981-01-01','2019-12-31'))
#Regional mean to create a time series
da = da.mean(dim=('lat','lon'))

#%%
##Step 3:  Sort the data from high to low values
da = da.sortby(da,ascending = False)
# Create a vector that convert the indices to values between 1 and #num_pts to values between 0 and 100
x = np.arange(1,len(da)+1)*100/len(da) #Adding 1 in np arange ensure the last value is 100% - important for small numbers of values
# Plot
plt.plot(x,da)
plt.xlabel("Chance of exceeding [%]")
plt.ylabel("Precipitation [mm/day]")
plt.show()