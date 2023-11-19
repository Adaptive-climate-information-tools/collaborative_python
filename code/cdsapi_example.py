#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 09:48:19 2023

@author: ellendyer
"""

# code example uses this dataset page: https://cds.climate.copernicus.eu/cdsapp#!/dataset/seasonal-monthly-single-levels?tab=form


import cdsapi

c = cdsapi.Client()

c.retrieve(
    'seasonal-monthly-single-levels',
    {
        'format': 'netcdf',
        'originating_centre': 'ecmwf',
        'system': '5',
        'variable': ['total_precipitation'],
        'product_type': [
            'hindcast_climate_mean', 'monthly_mean',
        ],
        'year': [
            '2010',
            '2011', '2012', '2013',
            '2014', '2015', '2016',
            '2017', '2018', '2019',
            '2020', '2021', '2022',
        ],
        'month': ['07','08','09'],
        'leadtime_month': '2',
        'area': [
            35, 20, -10,
            55,
        ],
    },
    'download.nc')