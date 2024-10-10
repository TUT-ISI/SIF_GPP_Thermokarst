"""
Script: This script downloads the required meteorological variables from the ECMWF ERA-5 products. 
The data is aggregated per year and month and stored in NetCDF.
Author: Lauri Hatakka
Version: 1.0
Date: Jun, 2023

"""
import cdsapi
import zipfile
import xarray as xr
import os
import shutil

SAVELOCATION = '/fmi/scratch/project_2005798/data/ERA5_Land_8day'
months = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']
# years = [] 2000 - 2020
years = list(range(2000,2021))
print('Starting the program!')
for year in years:
    print(f'###############################################')
    print(f'###############################################')
    print(f'STARTING YEAR {year}')
    print(f'###############################################')
    print(f'###############################################')   
    for month in months:

        c = cdsapi.Client()
        FILENAME = f'{SAVELOCATION}/tmp/{year}_{month}/{year}_{month}.netcdf.zip'
        if not os.path.isfile(FILENAME):
            print(f'The data for {year}/{month} was not found on disk, downloading!')
            os.makedirs(f'{SAVELOCATION}/tmp/{year}_{month}/')
            c.retrieve(
                'reanalysis-era5-land',
                {
                    'variable': [
                        'forecast_albedo', 'leaf_area_index_high_vegetation', 'leaf_area_index_low_vegetation',
                        'surface_latent_heat_flux', 'surface_net_solar_radiation', 'surface_pressure',
                        'surface_sensible_heat_flux', 'surface_solar_radiation_downwards', 'total_precipitation',
                    ],
                    'year': year,
                    'month': month,
                    'day': [
                        '01', '02', '03',
                        '04', '05', '06',
                        '07', '08', '09',
                        '10', '11', '12',
                        '13', '14', '15',
                        '16', '17', '18',
                        '19', '20', '21',
                        '22', '23', '24',
                        '25', '26', '27',
                        '28', '29', '30',
                        '31',
                    ],
                    'time': [
                        '00:00', '01:00', '02:00',
                        '03:00', '04:00', '05:00',
                        '06:00', '07:00', '08:00',
                        '09:00', '10:00', '11:00',
                        '12:00', '13:00', '14:00',
                        '15:00', '16:00', '17:00',
                        '18:00', '19:00', '20:00',
                        '21:00', '22:00', '23:00',
                    ],
                    'area': [
                        90, -180, 60,
                        180,
                    ],
                    'format': 'netcdf.zip',
                },
                FILENAME)
            print(f'Done!')

        else:
            print(f'The data for {year}/{month} was already downloaded!')

        tmplocation =f'{SAVELOCATION}/tmp/{year}_{month}' 

        if not os.path.isfile(f'{tmplocation}/data.nc'):
            print(f'The netcdf for {year}/{month} was not found, extracting from .zip!')
            with zipfile.ZipFile(FILENAME, 'r') as zip_ref:
                zip_ref.extractall(tmplocation) 
            print('Done!')
            
        else:
            print(f'Netcdf for {year}/{month} already existed!')
           
        
        print(f'Resampling {year}/{month} data to one day resolution! ')
        datafile = f'{SAVELOCATION}/tmp/{year}_{month}/data.nc'
        data = xr.open_dataset(datafile)  
        data = data.resample(time = '1D').mean(dim = 'time', skipna = True)
        data.to_netcdf(datafile)
        print(f'Done!')

        print(f'')
        print(f'###########')
        print(f'Aggregation for {year}/{month} is done and saved into a file data.nc in corresponding folder!')
        print(f'###########')
        print(f'')

    print(f'###############################################')
    print(f'All months for {year} done, starting aggregation')
    print(f'###############################################')


    data = xr.Dataset()
    
    for month in months:
        print(f'starting with month: {year}/{month}')
        datalocation = f'{SAVELOCATION}/tmp/{year}_{month}'
        data = xr.merge([data,xr.open_dataset(f'{datalocation}/data.nc', chunks = {'time' : 1, 'latitude' : 1000, 'longitude' : 1000})], join = 'outer')
        print(f'Merged!')
    
    print(f'Resampling year {year}!')
    data = data.resample(time = '8D').mean(dim = 'time', skipna = True)
    data.to_netcdf(f'{SAVELOCATION}/aggregateddata_{year}.nc')
    print(f'Done!')

    print(f'Removing files for {year}')
    shutil.rmtree(f'{SAVELOCATION}/tmp', ignore_errors = True)
    print(f'Done!')

print(f'EVERYTHING DONE!')
print(f'CLOSING THE PROGRAM :)')

    
