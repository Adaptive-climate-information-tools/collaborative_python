#%%
#Read in all the packages you need. 
#For this code all of the packages below are needed
import sys
import os
import glob
import xarray as xr
import pandas as pd
import numpy as np
from scipy import stats
from matplotlib import pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature
import matplotlib as mpl
mpl.rcParams['font.size'] = 10
mpl.rcParams['legend.fontsize'] = 10
mpl.rcParams['figure.titlesize'] = 10
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from matplotlib.ticker import PercentFormatter

#We will be taking a difference between future and historical time slices (slices
#of time that are the same length) so we will open each of these sets of model 
#runs separately:

#%%
##HISTORICAL CMIP6 READ IN
#Set your path to where the "cmip_files" folder is (we will learn to download these later)
path = '/Users/ellendyer/Desktop/Nov23_Python/collaborative_python/cmip_files/cmip_hist_nc_files/'
#All the files have different names but they start with the same string so you can use * to
#list them all
filename = 'ts_Amon_*'
#Try printing list_f to see what glob.glob does
list_f = glob.glob(os.path.join(path,filename))

#We will now loop through all of the model files to add them to one dataarray
#This is similar to how we loop through spreadsheet files too:
#We add them all to a list and then concatenate them - xarray does this with 
#the information about dimensions in the files
chl = []
daterange = pd.date_range(start='1980-01-01', end='2014-12-31', freq='MS')
#We pick a date range for our historical timeslice (we will use this date range
#to make sure all the model data has the same time stamps)
for m in list_f:
  #This will extract the name of the model from the file name
  mname = m.split("_")[-6]
  print(mname)
  #Convert from K to C - 273K
  inH = xr.open_dataset(m)['ts']-273
  #We don't really have to resample because we are using monthly data but
  #it is here in case you read in daily data
  inH = inH.resample(time="MS").mean()
  #We add a dimension for the model number because all the models
  #will be put in the same xarray and we will make the coordinate the
  #full name of the model
  inH = inH.expand_dims(dim="model")
  inH = inH.assign_coords(model=('model',[mname]))
  #Replace the time dimension with the date range above - now all the models
  #will have the same time stamps
  inH['time'] = daterange
  #Select a smaller regional section and interpolate to the same 
  #uniform grid for each model so math can be done with all the models
  inH = inH.sel(lat=slice(-10,26),lon=slice(21,61)).interp(lat=np.arange(-10,26,1),lon=np.arange(22,61,1), method="linear")
  #Save these nicely formatted individual model files to a new directory for use later
  inH.to_netcdf('../cmip_files/cmip_hist_nc_files/new/ts_Historical_'+mname+'.nc',mode='w')
  #Add the array for this model to the list of model arrays
  chl.append(inH)
#Concatenate all the invididual model arrays together into one arrah
cmipHin=xr.concat(chl,dim='model')

#print(cmipHin)

#%%
##SSP585 CMIP6 READ IN

#The routine to run in future model runs is basically the same as above except
#for the paths and date range

path = '/Users/ellendyer/Desktop/Nov23_Python/collaborative_python/cmip_files/cmip_sp585_nc_files/'
filename = 'ts_Amon_*'
list_f = glob.glob(os.path.join(path,filename))

chl = []
#daterange = pd.date_range(start='2015-01-01', end='2099-12-31', freq='MS')
daterange = pd.date_range(start='2066-01-01', end='2099-12-31', freq='MS')
for m in list_f:
    try:
        mname = m.split("_")[-6]
        print(mname)
        inH = xr.open_dataset(m)['ts']-273
        inH = inH.resample(time="MS").mean()
        inH = inH.sel(time=slice('2066-01-01','2099-12-31'))
        #inH = inH.sel(time=slice('2015-01-01','2099-12-31'))
        inH = inH.expand_dims(dim="model")
        inH = inH.assign_coords(model=('model',[mname]))
        inH['time'] = daterange
        inH = inH.sel(lat=slice(-10,26),lon=slice(21,61)).interp(lat=np.arange(-10,26,1),lon=np.arange(22,61,1), method="linear")
        inH.to_netcdf('../cmip_files/cmip_sp585_nc_files/new/ts_sp585_'+mname+'.nc',mode='w')
        chl.append(inH)
    except:
      print('model bad time')
cmipFin=xr.concat(chl,dim='model')


#print(cmipFin)

#%%

# models = list((set(cmipFin.model.values.tolist())).intersection(set(cmipHin.model.values.tolist())))
# cmipHin = cmipHin.sel(model=np.in1d(cmipHin['model'], models),drop=True)
# cmipFin = cmipFin.sel(model=np.in1d(cmipFin['model'], models),drop=True)

seas = 'OND'
seasT = 'October-December'
months = [10,11,12]

#select season
cmipH = cmipHin.sel(time=np.in1d(cmipHin['time.month'], months))


#select season
cmipF = cmipFin.sel(time=np.in1d(cmipFin['time.month'], months))


#Same models for historical as future
cmipH = cmipH.sel(model=np.in1d(cmipH['model'], list(cmipF.model.values)))


#%%
#Calculate some values with the historical and future model arrays

#ENSEMBLE MEAN

#Historical ensemble mean
cmipH_eM = cmipH.mean(dim='model')
#Historical ensemble std dev (model spread)
cmipH_eV = cmipH.std(dim='model').mean('time')

#Future ensemble mean
cmipF_eM = cmipF.mean(dim='model')
#Future ensemble std dev (model spread)
cmipF_eV = cmipF.std(dim='model').mean('time')

#Future change for each model
cmip_diff_models = cmipF.mean('time') - cmipH.mean('time')
#Ensemble mean future change
cmip_diff_ens = (cmipF.mean('time') - cmipH.mean('time')).mean('model')


#Get the agreement in the sign change among models to 
#plot as a hash on top of the contour plots of rainfall values

#Sign of future change for each model (either 1 or -1)
cmip_diff_sign = cmip_diff_models/np.abs(cmip_diff_models)

#For models where the sign of their future change isn't 
#the same as the sign of the ensemble average change set value to zero
cmip_diff_sign_models = cmip_diff_sign.where(np.sign(cmip_diff_sign) == np.sign(cmip_diff_ens),0)

#Sum up the model field so that there is a count of models that agree
cmip_diff_sign_agree = cmip_diff_sign_models.sum('model')
#Total number of possible models in ensemble
num_models = len(cmipF.model)
#Percent of models that agree with the ensemble sign at each point
cmip_diff_sign_agree_percent = cmip_diff_sign_agree/num_models
#Drop points where 80% of models don't agree on the sign of ensemble change
thresh_agree = cmip_diff_sign_agree_percent.where(np.abs(cmip_diff_sign_agree_percent) > 0.8,drop=True)
#Turn this into 1s and nans for stippling plots
thresh_agree_point = thresh_agree/thresh_agree

#%%

#Plot the difference in the ensemble mean and overlay with hash lines
#where %80 of models agree on the sign change

#Set up the plot axes and features with cartopy
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()
ax.add_feature(cartopy.feature.BORDERS,lw=2)
ax.add_feature(cartopy.feature.COASTLINE,lw=2)
ax.add_feature(cartopy.feature.LAKES,lw=1)
ax.add_feature(cartopy.feature.RIVERS,lw=1)

#Select a python colormap (https://matplotlib.org/stable/gallery/color/colormap_reference.html)
#and we also modify it to set values over the max value to be a contrast colour
cmap = plt.get_cmap('YlOrRd').copy()
cmap.set_extremes(over='pink')

#Plot the ensemble difference using a colourmesh - you could use 
#a filled contour plot too
cmip_diff_ens.plot.pcolormesh(ax=ax,cmap=cmap,transform=ccrs.PlateCarree()
               ,vmin=0,vmax=5
               ,extend='both'
               ,robust=True,cbar_kwargs={'orientation':'horizontal'
                                         ,'fraction':0.04
                                         ,'pad':0.09
                                         ,'label': ""})
#Plot the sign agreement with a filled contour where the fill is a 
#hatch patter "////"
thresh_agree_point.plot.contourf(levels=[0.5,1.5],hatches=["////"],
                                  alpha=0.0,
                                  colors='none',add_colorbar=False)

#Set the plot title - using the season name set earlier in the code
ax.set_title(seasT+' temperature change',fontsize=9)

#Set the grid lines on the plot (style and location)
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='black', alpha=0.5, linestyle='--')
gl.top_labels = False
gl.right_labels = False
gl.xlocator = mticker.FixedLocator([25,35,45,55])
gl.ylocator = mticker.FixedLocator([-10,0,10,20])
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER

#Set the spatial extent of the plot [West, East, South, North]
ax.set_extent([22, 60, -8, 20])

#Annotate the colourbar with a unit (x,y,label,textsize)
plt.text(22,-13, '$^\circ$C', dict(size=11))

#Save the figure - set your own path!
#This uses the season acronym set earlier in the code
plt.savefig('../plots/cmip_Fdiff_ts_'+seas+'.png',bbox_inches='tight',dpi=300)

#Preview the plot
plt.show()
plt.clf()


#%%

#Use the same plotting routine as above but plot the spread (std dev) of the 
#distribution of model futures

ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()
ax.add_feature(cartopy.feature.BORDERS,lw=2)
ax.add_feature(cartopy.feature.COASTLINE,lw=2)
ax.add_feature(cartopy.feature.LAKES,lw=1)
ax.add_feature(cartopy.feature.RIVERS,lw=1)
cmap = plt.get_cmap('magma_r').copy()
cmap.set_extremes(over='darkgrey')

cmipF_eV.plot.pcolormesh(ax=ax,cmap=cmap,transform=ccrs.PlateCarree()
               ,vmin=0,vmax=5
               ,extend='both'
               ,robust=True,cbar_kwargs={'orientation':'horizontal'
                                         ,'fraction':0.04
                                         ,'pad':0.09
                                         ,'label': ""})

ax.set_title(seasT+'\n spread in future temperature',fontsize=9)

gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='black', alpha=0.5, linestyle='--')
gl.top_labels = False
gl.right_labels = False
gl.xlocator = mticker.FixedLocator([25,35,45,55])
gl.ylocator = mticker.FixedLocator([-10,0,10,20])
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
ax.set_extent([22, 60, -8, 20])
plt.text(22,-13, '$^\circ$C', dict(size=11))
plt.savefig('cmip_F_ts_spread_'+seas+'.png',bbox_inches='tight',dpi=300)
plt.show()
plt.clf()



#%%

#Use the same plotting routine as above but plot the spread (std dev) of the 
#distribution of model futures

ax = plt.axes()
cmipH.plot.hist(bins=60,density=True,label='Historical 1980-2015')
cmipF.plot.hist(bins=60,alpha=0.5,density=True,label='SP585 2066-2100')
ax.legend()
ax.yaxis.set_major_formatter(PercentFormatter(xmax=1))
ax.set_xlabel('monthly temperature $^\circ$C')
ax.set_ylabel('percent of total events in the period')
plt.savefig('../plots/cmip_ts_dist_'+seas+'.png',bbox_inches='tight',dpi=300)
plt.show()
plt.clf()


#%%

#Use the same plotting routine as above but plot timeseries of SP585 temperatures
#to the end of the century

fig, ax = plt.subplots(figsize=(10, 5))
for m in cmipF.model.values:
    print(m)
    seas_tl = cmipF.sel(model=m).mean(dim=('lat','lon'))
    seas_tl = seas_tl.groupby('time.year').mean('time')
    trend = stats.linregress(seas_tl.year.values, seas_tl.values)
    slope = trend[0]
    seas_tl.plot(label=m+' ('+'{:.2f}'.format(slope)+'$^\circ$C/year)')
ax.legend(fontsize=8)
ax.set_xlabel('Year')
ax.set_ylabel('Temperature $^\circ$C')
ax.set_title('Temperature trend')
plt.savefig('../plots/cmip_ts_tseries_'+seas+'.png',bbox_inches='tight',dpi=300)
plt.show()
plt.clf()



