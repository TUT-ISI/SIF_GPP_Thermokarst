"""
Script: This code reads the data (SIF and GPP) from the original source and temporal resolution and store it in a NetCDF on a 8-day interval. 
# Also, it computes the slope between the SIF and the different GPP products. Everything is resampled to the same spatial and temporal 
# resolution.
Author: Lauri Hatakka
Version: 1.0
Date: Jun, 2023

"""
import os
import code
import sys
import numpy as np
from netCDF4 import Dataset, date2num
from glob import glob
from datetime import datetime, timedelta
import xarray as xr
import rasterio as rio
import zipfile


def keyboard(banner=None):
    ''' Function that mimics the matlab keyboard command '''
    # use exception trick to pick up the current frame
    try:
        raise None
    except Exception:
        frame = sys.exc_info()[2].tb_frame.f_back
    # evaluate commands in current namespace
    namespace = frame.f_globals.copy()
    namespace.update(frame.f_locals)
    try:
        code.interact(banner=banner, local=namespace)
    except SystemExit:
        return

# s.sys.argv catches the argument when calling the python script 
# in command line, i.e. the inputof a function
year = int(os.sys.argv[1])
print('Processing year {year:d}'.format(year=year))

local_path = '../data/SIFGPPX/'
outputfilename = 'SIFGPPx_{year:d}.nc4'.format(year=year)
outputfilename = os.path.join(local_path,outputfilename)

os.makedirs(local_path, exist_ok=True)
#outputfilename = 'data/SIFGPPX/SIFGPPxSlope_{year:d}.nc4'.format(year=year)
if os.path.isfile(outputfilename):
    print('Done already!')
    sys.exit(0)

# From the GOSIF database
GOSIFfiles = glob('/fmi/scratch/project_2005798/data/GOSIF_v2_8day/*.tif.gz')
GOSIFdates = sorted([int(x[-14:-7]) for x in GOSIFfiles])

try:
    datasetGPPFluxSat = xr.open_mfdataset('/fmi/scratch/project_2005798/data/GPP_FluxSat_daily_v2/GPP_FluxSat_daily_v2_{year:d}*.nc4'.format(year=year), combine='by_coords')
except Exception:
    print('ERROR IN LOADING GPP')
    sys.exit(0)

lat, lon = None, None
GPP, SIF = [], []
times = []
for ii, GOSIFdaterange0 in enumerate(GOSIFdates):
    N0 = datetime(int(str(GOSIFdaterange0)[:4]), 1, 1) + timedelta(days=int(str(GOSIFdaterange0)[-3:]) - 1)
    if N0.year != year:
        continue

    print(N0)

    if ii < len(GOSIFdates) - 1:
        this_date = GOSIFdates[ii + 1] - 1
        N1 = datetime(int(str(this_date)[:4]), 1, 1) + timedelta(days=int(str(this_date)[-3:]))
    else:
        N1 = N0 + timedelta(days=8)
        if N1.year != N0.year:
            N1 = datetime(N1.year, 1, 1)
    print('    ', N0, N1)
    times.append(N0)

    GOSIFfilename = '/vsigzip//fmi/scratch/project_2005798/data/GOSIF_v2_8day/GOSIF_{}.tif.gz'.format(GOSIFdaterange0)
    #GOSIFfilename = '/vsigzip/data/GOSIF_v2_8day/GOSIF_{}.tif.gz'.format(GOSIFdaterange0)
    GOSIFdata = np.array(rio.open(GOSIFfilename).read()).astype(float)[0]
    GOSIFdata[GOSIFdata == 32767] = np.nan  # water bodies
    GOSIFdata[GOSIFdata == 32766] = np.nan  # lands under snow/ice throughout the year
    SIF.append(GOSIFdata[::-1, :] * 0.0001)

    GPPdata = datasetGPPFluxSat.sel(time=slice(N0, N1)).mean('time')
    GPP.append(GPPdata['GPP'].values)

    if lat is None:
        lat = GPPdata.lat.values
        lon = GPPdata.lon.values

datasetGPPFluxSat.close()

SIF = np.stack(SIF, axis=0)
GPP = np.stack(GPP, axis=0)



# From the Sun Shade database
local_path = '/fmi/scratch/project_2005798/data/SIF_Sun_Shade/'
GPPfiles = glob(''.join([local_path, '{year:d}/8day/GPP_v21_{year:d}_*.tif'.format(year=year)]))
GPPdates = sorted([int(x[-7:-4]) for x in GPPfiles])

print(''.join([local_path, '{year:d}.zip'.format(year=year)]))
with zipfile.ZipFile(''.join([local_path, '{year:d}.zip'.format(year=year)]), 'r') as zip_ref:
    zip_ref.extractall(local_path)   
#try:
 #   with zipfile.ZipFile(''.join([local_path, '{year:d}.zip'.format(year=year)]), 'r') as zip_ref:
  #       zip_ref.extractall(local_path)    
#except Exception:
#    print('ERROR IN LOADING GPP_SUN_SHADE data')
#    sys.exit(0)
#    a = oooooooo


latG, lonG = None, None
GPPSun, GPPShade, GPPtotal = [], [], []
times = []
scale_factor = 0.01 # Information from https://doi.org/10.1038/s41597-022-01309-2
# Note that the scale factor for monthly data is different: 0.1
scale_day = 1/8 # GPP dataset is provided by 8 day-1. So, to make it daily it should be 1/8
   
for ii, GPPdaterange0 in enumerate(GPPdates):
    N0 = datetime(year,1,1) + timedelta(days=int(str(GPPdaterange0)[-3:]) - 1)

    if N0.year != year:
        continue
    
    if ii < len(GPPdates) - 1:
        this_date = GPPdates[ii + 1] - 1
        N1 = datetime(year,1,1) + timedelta(days=int(str(this_date)[-3:]))
    else:
        N1 = N0 + timedelta(days=8)
        if N1.year != N0.year:
            N1 = datetime(N1.year, 1, 1)
    print('    ', N0, N1)
    times.append(N0)
    
    GPPShadefilename =''.join([local_path, '{year:d}/8day/Shade_GPP_v21_{year:d}_{doy:03d}.tif'.format(year=year, doy=GPPdaterange0)])
    GPPSunfilename =''.join([local_path, '{year:d}/8day/Sun_GPP_v21_{year:d}_{doy:03d}.tif'.format(year=year, doy=GPPdaterange0)])
    GPPtotalfilename =''.join([local_path, '{year:d}/8day/GPP_v21_{year:d}_{doy:03d}.tif'.format(year=year, doy=GPPdaterange0)])
    print(GPPSunfilename)


    if lat is None:        
        dataset = rio.open(GPPShadefilename)
        interval = abs((dataset.bounds.left-dataset.bounds.right)/dataset.shape[1])
        lonG                   = np.linspace(dataset.bounds.left+(interval/2), dataset.bounds.right+(interval/2), dataset.shape[1])
        latG                   = np.linspace(dataset.bounds.bottom+(interval/2), dataset.bounds.top+(interval/2), dataset.shape[0])
        #[lon,lat] = np.meshgrid(x,y)
       
    
    GPPShadedata = np.array(rio.open(GPPShadefilename).read()).astype(float)[0]
    # keyboard()
    GPPShadedata = GPPShadedata*scale_factor*scale_day
    #GPPShadedata[GPPShadedata == 0] = np.nan  # water bodies
    GPPShade.append(GPPShadedata[::-1, :])
    
    GPPSundata = np.array(rio.open(GPPSunfilename).read()).astype(float)[0]
    GPPSundata = GPPSundata*scale_factor*scale_day
    #GPPSundata[GPPSundata == 0] = np.nan  # water bodies
    GPPSun.append(GPPSundata[::-1, :])
    
    GPPtotaldata = np.array(rio.open(GPPtotalfilename).read()).astype(float)[0]
    GPPtotaldata = GPPtotaldata*scale_factor*scale_day
    #GPPtotaldata[GPPtotaldata == 0] = np.nan  # water bodies
    GPPtotal.append(GPPtotaldata[::-1, :])


GPPSun   = np.stack(GPPSun, axis=0)
GPPShade = np.stack(GPPShade, axis=0)
GPPtotal = np.stack(GPPtotal, axis=0)


# Check if the latG and lat and lonG and lon are the same 

#if latG==lat:
 #   print('LAT is equal in the two datasets')
#if lonG==lon:
#    print('LON is equal in the two datasets')


# -------------------------------------------------------
# -------------------------------------------------------
# write the netCDF file
# -------------------------------------------------------
# -------------------------------------------------------
ncout = Dataset(outputfilename, 'w', format='NETCDF4')

# Add some attributes
ncout.History = 'File generated on {} (UTC) by {}'.format(datetime.utcnow().strftime('%c'), os.path.basename(__file__))
ncout.Description = '8 day combined GOSIF_v2_8day and GPP_FluxSat_daily_v2, and the  GPPsun, GPPShade, and GPPtotal. Dataset source https://doi.org/10.1038/s41597-022-01309-2'

# create dimensions
ncout.createDimension('time', len(times))
ncout.createDimension('latitude', len(lat))
ncout.createDimension('longitude', len(lon))

# save coordinates
ncout_latitude = ncout.createVariable('latitude', 'f4', ('latitude',))
ncout_latitude[:] = lat
ncout_longitude = ncout.createVariable('longitude', 'f4', ('longitude',))
ncout_longitude[:] = lon

# create & save time
units = 'seconds since 1980-01-01'
ncout_time = ncout.createVariable('time', 'i4', ('time',))
ncout_time[:] = date2num(times, units)
ncout_time.units = units

# save data
ncout_dist = ncout.createVariable('GPP', 'f4', ('time', 'latitude', 'longitude'))
ncout_dist.long_name = 'Gross primary productivity of biomass expressed as carbon'
ncout_dist.units = 'g m-2 d-1'
ncout_dist[:] = GPP

ncout_dist = ncout.createVariable('SIF', 'f4', ('time', 'latitude', 'longitude'))
ncout_dist.long_name = 'Solar induced fluorescence'
ncout_dist.units = 'W m−2 um−1 sr−1'
ncout_dist[:] = SIF

ncout_dist = ncout.createVariable('GPPSun', 'f4', ('time', 'latitude', 'longitude'))
ncout_dist.long_name = 'Direct sun induced Gross primary productivity of biomass expressed as carbon '
ncout_dist.units = 'g C m-2 day-1'
ncout_dist[:] = GPPSun

ncout_dist = ncout.createVariable('GPPShade', 'f4', ('time', 'latitude', 'longitude'))
ncout_dist.long_name = 'Direct sun induced Gross primary productivity of biomass expressed as carbon '
ncout_dist.units = 'g C m-2 day-1'
ncout_dist[:] = GPPShade

ncout_dist = ncout.createVariable('GPPtotal', 'f4', ('time', 'latitude', 'longitude'))
ncout_dist.long_name = 'Direct sun induced Gross primary productivity of biomass expressed as carbon '
ncout_dist.units = 'g C m-2 day-1'
ncout_dist[:] = GPPtotal

ncout.close()
