"""
Script: Read from the ARCLIM database the 2mt, convert into C and crop the latitude region of interest
Author: Lauri Hatakka
Version: 1.0
Date: Jun, 2023

Data available in 
https://www.nature.com/articles/s41597-023-01959-w
For daily resolution data ask the authors.
"""
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt



years = list(range(2000,2021))
months = list(range(1,13))
TEST = False
if TEST:
    years = [2000]

for year in years:
    print(f'Starting {year}')
    dataset = xr.Dataset()
    outpath = f'/yourpath/ARCLIM/self_calculated/2t_DMEA_{year}.nc'
    for month in months:
        print(f'Starging {month}')
        da_t2mean = xr.open_dataset(f'/ARCLIM_path/daily_from_allas/2m_temperature_DMEA_era5Land_{year}{month:02}.nc')

        # Convert temperatures to celsius
        da_t2mean['t2m'] = da_t2mean.t2m - 273.15 

        # Drop uninteresting latitudes
        ds_2t_month = da_t2mean.sel(latitude = slice(90.0, 60.0))    
        
        dataset = xr.merge([dataset, ds_2t_month])

    # Aggregate yearly data into 8-day interval, which will be used in the future analysis
    dataset = dataset.resample(time = '8D').mean(dim='time', skipna = True)
    print(dataset)
       

    # save the data as a netcdf file
    dataset.to_netcdf(outpath, format='NETCDF4')
    print(f'######################')
    print(f'{year} done!')
    print(f'######################')