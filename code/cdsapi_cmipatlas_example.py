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
    'projections-climate-atlas',
    {
        'format': 'zip',
        'origin': 'CMIP6',
        'experiment': 'historical',
        'domain': 'global',
        'period': '1850-2014',
        'variable': 'monthly_mean_of_daily_accumulated_precipitation',
    },
    'cmipatlas_download.zip')