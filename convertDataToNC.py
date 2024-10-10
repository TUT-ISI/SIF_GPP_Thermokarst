"""
Script: Read and store the GPP data from Fluxsat_v2 database and SIF from GOSIF databases in a NetCDF output format.
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


year = int(os.sys.argv[1])
print('Processing year {year:d}'.format(year=year))
os.makedirs('data/nc', exist_ok=True)
outputfilename = 'data/nc/SIFGPP_{year:d}.nc4'.format(year=year)
if os.path.isfile(outputfilename):
    print('Done already!')
    sys.exit(0)

GOSIFfiles = glob('data/GOSIF_v2_8day/*.tif.gz')
GOSIFdates = sorted([int(x[-14:-7]) for x in GOSIFfiles])

try:
    datasetGPPFluxSat = xr.open_mfdataset('data/GPP_FluxSat_daily_v2/GPP_FluxSat_daily_v2_{year:d}*.nc4'.format(year=year), combine='by_coords')
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

    GOSIFfilename = '/vsigzip/data/GOSIF_v2_8day/GOSIF_{}.tif.gz'.format(GOSIFdaterange0)
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

GPP = np.stack(GPP, axis=0)
SIF = np.stack(SIF, axis=0)

# write the netCDF file
ncout = Dataset(outputfilename, 'w', format='NETCDF4')

# Add some attributes
ncout.History = 'File generated on {} (UTC) by {}'.format(datetime.utcnow().strftime('%c'), os.path.basename(__file__))
ncout.Description = '8 day combined GOSIF_v2_8day and GPP_FluxSat_daily_v2'

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

ncout.close()
