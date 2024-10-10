"""
Script: This script prepares and merge all the inpyt data required to perform the SHAP analysis for GPP. 
Author: Lauri Hatakka
Version: 1.0
Date: Jun, 2023

"""
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import datetime
import os
import glob
import pandas as pd
import dask



FILELOCATION = '/fmi/scratch/project_2005798/data/'
years = list(range(2000,2021))
month_groups = [[4,5],[6,7],[8,9]]


for months_wanted in month_groups:
    for year in years:
        outputdir = f'{FILELOCATION}PreprocessedData/ERAGPP/{months_wanted}/'
        os.makedirs(outputdir, exist_ok = True)
        try:
            ERA_data = xr.open_dataset(os.path.join(FILELOCATION, f'ERA5_Land_8day/aggregateddata_{year}.nc'), chunks = 'auto')
            ERA_data = ERA_data.isel(time = ERA_data.time.dt.month.isin(months_wanted))
            VPD_data = xr.open_dataset(os.path.join(FILELOCATION, f'ARCLIM/self_calculated/vpd_DMEA_{year}.nc'), chunks = 'auto')
            VPD_data = VPD_data.isel(time=VPD_data.time.dt.month.isin(months_wanted))
            TEMP_data = xr.open_dataset(os.path.join(FILELOCATION, f'ARCLIM/self_calculated/2t_DMEA_{year}.nc'), chunks = 'auto')
            TEMP_data = TEMP_data.isel(time = TEMP_data.time.dt.month.isin(months_wanted))

            SIFGPP_data = xr.open_dataset(os.path.join(FILELOCATION, f'GPPSunShadeNC_8day/SIFGPPSunShade_{year}.nc4'), chunks = 'auto')
            SIFGPP_data = SIFGPP_data.isel(time=SIFGPP_data.time.dt.month.isin(months_wanted))
            SIFGPP_data =SIFGPP_data.sel(latitude = slice(60,91))

        except FileNotFoundError as e:
            print(e)
            print(f'All the files for year {year} were not found, skipping!')
            continue
        ERA_data = ERA_data.to_dataframe().dropna(how='all').to_xarray()
        VPD_data = VPD_data.to_dataframe().dropna(how='all').to_xarray()
        TEMP_data = TEMP_data.to_dataframe().dropna(how = 'all').to_xarray()
        SIFGPP_data = SIFGPP_data.coarsen(latitude=2, longitude=2 , coord_func = 'min', boundary = 'trim').mean()
        SIGFPP_data = SIFGPP_data.to_dataframe().dropna(how='any').to_xarray()
        SIFGPP_data = SIFGPP_data.reindex_like(ERA_data, method = 'nearest')
        
        dataset = xr.merge([ERA_data, SIFGPP_data, VPD_data, TEMP_data], join='inner')
        dataset.to_netcdf(f'{outputdir}{year}.nc')
        print(dataset)