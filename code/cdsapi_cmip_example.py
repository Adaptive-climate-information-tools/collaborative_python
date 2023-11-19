#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 09:48:19 2023

@author: ellendyer
"""

# code example uses this dataset page: https://cds.climate.copernicus.eu/cdsapp#!/dataset/projections-cmip6?tab=form


import cdsapi
import zipfile

for var in ['precipitation','surface_temperature']:
    for mod in ['hadgem3_gc31_mm','cesm2', 'canesm5', 'cnrm_cm6_1_hr','mpi_esm1_2_lr',
                'mri_esm2_0','ukesm1_0_ll','bcc_csm2_mr','noresm2_mm','cmcc_cm2_sr5']:
        c = cdsapi.Client()
        
        c.retrieve(
            'projections-cmip6',
            {
                'format': 'zip',
                'temporal_resolution': 'monthly',
                'experiment': 'historical',
                'variable': var,
                'model': mod,
                'month': [
                    '01', '02', '03',
                    '04', '05', '06',
                    '07', '08', '09',
                    '10', '11', '12',
                ],
                'year': [
                    '1980', '1981', '1982',
                    '1983', '1984', '1985',
                    '1986', '1987', '1988',
                    '1989', '1990', '1991',
                    '1992', '1993', '1994',
                    '1995', '1996', '1997',
                    '1998', '1999', '2000',
                    '2001', '2002', '2003',
                    '2004', '2005', '2006',
                    '2007', '2008', '2009',
                    '2010', '2011', '2012',
                    '2013', '2014',
                ],
                'area': [
                    40, 12, -10,
                    70,
                ],
            },
            'temp.zip')
        with zipfile.ZipFile('temp.zip', 'r') as zip_ref:
            for i in zip_ref.namelist():
                if i[-3:]=='.nc':
                    zip_ref.extract(i,path='cmip_hist_nc_files')        
        
for var in ['precipitation','surface_temperature']:
    for mod in ['hadgem3_gc31_mm','cesm2', 'canesm5', 'cnrm_cm6_1_hr','mpi_esm1_2_lr',
                'mri_esm2_0','ukesm1_0_ll','bcc_csm2_mr','noresm2_mm','cmcc_cm2_sr5']:
        c = cdsapi.Client()
        
        c.retrieve(
            'projections-cmip6',
            {
                'format': 'zip',
                'temporal_resolution': 'monthly',
                'experiment': 'ssp5_8_5',
                'variable': var,
                'model': mod,
                'month': [
                    '01', '02', '03',
                    '04', '05', '06',
                    '07', '08', '09',
                    '10', '11', '12',
                ],
                'year': [
                    '2015', '2016', '2017',
                    '2018', '2019', '2020',
                    '2021', '2022', '2023',
                    '2024', '2025', '2026',
                    '2027', '2028', '2029',
                    '2030', '2031', '2032',
                    '2033', '2034', '2035',
                    '2036', '2037', '2038',
                    '2039', '2040', '2041',
                    '2042', '2043', '2044',
                    '2045', '2046', '2047',
                    '2048', '2049', '2050',
                    '2051', '2052', '2053',
                    '2054', '2055', '2056',
                    '2057', '2058', '2059',
                    '2060', '2061', '2062',
                    '2063', '2064', '2065',
                    '2066', '2067', '2068',
                    '2069', '2070', '2071',
                    '2072', '2073', '2074',
                    '2075', '2076', '2077',
                    '2078', '2079', '2080',
                    '2081', '2082', '2083',
                    '2084', '2085', '2086',
                    '2087', '2088', '2089',
                    '2090', '2091', '2092',
                    '2093', '2094', '2095',
                    '2096', '2097', '2098',
                    '2099',
                ],
                'area': [
                    40, 12, -10,
                    70,
                ],
            },
            'temp.zip')
        with zipfile.ZipFile('temp.zip', 'r') as zip_ref:
            for i in zip_ref.namelist():
                if i[-3:]=='.nc':
                    zip_ref.extract(i,path='cmip_sp585_nc_files')        

        
        
        
     