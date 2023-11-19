#%%
import sys
import os
import glob
import xarray as xr
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature
import geopandas as gpd
from shapely.geometry import mapping
import datetime
import seaborn as sns
from scipy.stats import spearmanr
import matplotlib as mpl
mpl.rcParams['font.size'] = 10
mpl.rcParams['legend.fontsize'] = 10
mpl.rcParams['figure.titlesize'] = 10
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER


#%%
##HISTORICAL CMIP6 READ IN
path = '/Volumes/blue_wd/cmip6_museum/'
filename = 'pr_mods_Historical_*.nc'
list_f = glob.glob(os.path.join(path,filename))

chl = []
daterange = pd.date_range(start='1980-01-01', end='2014-12-31', freq='MS')
#cmipH = xr.Dataset()
for m in list_f:
  #print(m)
  mname=m[49:-3]
  print(mname)
  inH = xr.open_dataset(m)['pr']
  inH = inH.resample(time="MS").mean()
  inH = inH.expand_dims(dim="model")
  inH = inH.assign_coords(M=('model',[mname]))
  inH['time'] = daterange
  inH.sel(lat=slice(-12,22),lon=slice(22,53)).to_netcdf('/Volumes/blue_wd/meron_files/pr_Historical_'+mname+'.nc')
  chl.append(inH)
  #print(inH.time)
cmipHin=xr.concat(chl,dim='model')

print(cmipHin)

#%%
##SSP585 CMIP6 READ IN
path = '/Volumes/blue_wd/cmip6_museum/'
filename = 'pr_mods_SSP585_*.nc'
list_f = glob.glob(os.path.join(path,filename))

chl = []
daterange = pd.date_range(start='2066-01-01', end='2100-12-31', freq='MS')
#cmipH = xr.Dataset()
for m in list_f:
  try:
      #print(m)
      mname=m[45:-3]
      print(mname)
      inH = xr.open_dataset(m)['pr']
      inH = inH.resample(time="MS").mean()
      inH = inH.expand_dims(dim="model")
      inH = inH.assign_coords(M=('model',[mname]))
      #print(len(inH.time))
      inH['time'] = daterange
      inH.sel(lat=slice(-12,22),lon=slice(22,53)).to_netcdf('/Volumes/blue_wd/meron_files/pr_SSP585_'+mname+'.nc')
      chl.append(inH)
  except:
      print('model bad time')
cmipFin=xr.concat(chl,dim='model')


print(cmipFin)

#%%

seas = 'OND'
seasT = 'October-December'
months = [10,11,12]

#select season
cmipH = cmipHin.sel(time=np.in1d(cmipHin['time.month'], months))


#select season
cmipF = cmipFin.sel(time=np.in1d(cmipFin['time.month'], months))


#Same models for historical as future
cmipH = cmipH.sel(model=np.in1d(cmipH['M'], list(cmipF.M.values)))


#%%
#ENSEMBLE MEAN
#Historical ensemble mean
cmipH_eM = cmipH.mean(dim='model')
#Historical ensemble std dev (model spread)
cmipH_eV = cmipH.std(dim='model')

#Future ensemble mean
cmipF_eM = cmipF.mean(dim='model')
#Future ensemble std dev (model spread)
cmipF_eV = cmipF.std(dim='model')

#Future change for each model
cmip_diff_models = cmipF.mean('time') - cmipH.mean('time')
#Ensemble mean future change
cmip_diff_ens = (cmipF.mean('time') - cmipH.mean('time')).mean('model')

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
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()
ax.add_feature(cartopy.feature.BORDERS,lw=2)
ax.add_feature(cartopy.feature.COASTLINE,lw=2)
ax.add_feature(cartopy.feature.LAKES,lw=1)
ax.add_feature(cartopy.feature.RIVERS,lw=1)
#cmap = plt.get_cmap('Greens').copy()
#cmap.set_extremes(under='brown')
cmip_diff_ens.plot(ax=ax,cmap=plt.cm.RdYlBu,transform=ccrs.PlateCarree()
               ,vmin=-1,vmax=1
               ,extend='both'
               ,robust=True,cbar_kwargs={'orientation':'horizontal'
                                         ,'fraction':0.04
                                         ,'pad':0.09
                                         ,'label': ""})
thresh_agree_point.plot.contourf(levels=[0.5,1.5],hatches=["////"],
                                 alpha=0.0,
                                 colors='none',add_colorbar=False)
ax.set_title(seasT+' rainfall change',fontsize=9)


gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='black', alpha=0.5, linestyle='--')
gl.top_labels = False
gl.right_labels = False
gl.xlocator = mticker.FixedLocator([25,35,45,55])
gl.ylocator = mticker.FixedLocator([-10,0,10,20])
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
ax.set_extent([22, 60, -15, 22])
plt.text(15,-21, 'mm/day', dict(size=9))
plt.savefig('cmip_diff_'+seas+'.png',bbox_inches='tight',dpi=300)
plt.show()
plt.clf()

