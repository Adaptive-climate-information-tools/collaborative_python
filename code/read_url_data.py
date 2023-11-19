import sys
import xarray as xr
import netCDF4
import pandas as pd
#------------------
import requests

def fix_calendar(ds, timevar):
    if ds[timevar].attrs['calendar'] == '360':
        ds[timevar].attrs['calendar'] = '360_day'
    return ds

# Download NMME model forecast or hindcast runs
url = "https://iridl.ldeo.columbia.edu/SOURCES/.Models/.NMME/.NCAR-CESM1/.FORECAST/.MONTHLY/.prec/dods"

ds = xr.open_dataset(url,decode_times=False)
print(ds)
ds = fix_calendar(ds,'S')
ds = xr.decode_cf(ds)
ds = ds.sel(S=ds.S.dt.month.isin([3,4,5]))
print(ds)
ds = ds.rename({'S':'forecast_reference_time'})
ds = ds.rename({'X':'lon'})
ds = ds.rename({'Y':'lat'})
ds = ds.rename({'M':'realization'})
ds = ds.rename({'L':'forecast_period'})
ds = ds.sel(lat=slice(4,15),lon=slice(32,50))
print(ds)


# Download FEWSNET dekadal data
url = "https://iridl.ldeo.columbia.edu/SOURCES/.USGS/.EROS/.FEWS/.dekadal/.EAF/.Belg/.dw/dods"

ds = xr.open_dataset(url,decode_times=True)
print(ds)
ds = ds.sel(T=ds.T.dt.month.isin([3,4,5]))
print(ds)
ds = ds.rename({'T':'time'})
ds = ds.rename({'X':'lon'})
ds = ds.rename({'Y':'lat'})
ds = ds.sel(lat=slice(4,15),lon=slice(32,50))
print(ds)


# Download CHIRPS rainfall monthly
url = "https://iridl.ldeo.columbia.edu/SOURCES/.UCSB/.CHIRPS/.v2p0/.monthly/.global/.precipitation/dods"

ds = xr.open_dataset(url,decode_times=False)
print(ds)
ds = fix_calendar(ds,'T')
ds = xr.decode_cf(ds)
ds = ds.sel(T=ds.T.dt.month.isin([3,4,5]))
print(ds)
ds = ds.rename({'T':'time'})
ds = ds.rename({'X':'lon'})
ds = ds.rename({'Y':'lat'})
ds = ds.sel(lat=slice(15,4),lon=slice(32,50))
print(ds)

# Download CHIRTS temperature daily (also use this method for daily CHIRPS)
url = "https://iridl.ldeo.columbia.edu/SOURCES/.UCSB/.CHIRTS/.v1.0/.daily/.global/.0p05/dods"

ds = xr.open_dataset(url)
print(ds.T)
pdtime = pd.to_datetime(ds.T,origin='julian',unit='D')
ds['T'] = pdtime
print(ds)
ds = ds.sel(T=ds.T.dt.month.isin([3,4,5]))
print(ds)
ds = ds.rename({'T':'time'})
ds = ds.rename({'X':'lon'})
ds = ds.rename({'Y':'lat'})
ds = ds.sel(lat=slice(4,15),lon=slice(32,50))
print(ds)

# Download TAMSAT rainfall daily
url = "https://iridl.ldeo.columbia.edu/SOURCES/.Reading/.Meteorology/.TAMSAT/.v3p1/.daily/.rfe/dods"

ds = xr.open_dataset(url)
print(ds.T)
pdtime = pd.to_datetime(ds.T,origin='julian',unit='D')
ds['T'] = pdtime
print(ds)
ds = ds.sel(T=ds.T.dt.month.isin([3,4,5]))
print(ds)
ds = ds.rename({'T':'time'})
ds = ds.rename({'X':'lon'})
ds = ds.rename({'Y':'lat'})
ds = ds.sel(lat=slice(4,15),lon=slice(32,50))
print(ds)

# # Download TAMSAT rainfall daily
# url = "https://iridl.ldeo.columbia.edu/SOURCES/.ECMWF/.S2S/.NCEP/.reforecast/.perturbed/.sfc_precip/.tp/dods"

# session = requests.Session()
# session.auth = ('ellen.dyer@ouce.ox.ac.uk','SillyCat')
# store = xr.backends.PydapDataStore.open(url,session=session)
# ds = xr.open_dataset(store)
# print(ds)

# from pydap.client import open_url
# from pydap.cas.urs import setup_session
# session = setup_session(username='ellen.dyer@ouce.ox.ac.uk',password='SillyCat',check_url='https://iridl.ldeo.columbia.edu/auth/login')
# ds = open_url(url,session)

