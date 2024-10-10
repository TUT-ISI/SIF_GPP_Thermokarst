"""
Script Name: convertPeatlandTifToNC.py
Description: This script convert information from the Histel and Histosol datasets 
             from .tiff to .nc The original data can be found as supplemental data
             for the scientific paper Hugelius et al 2020 (https://www.pnas.org/cgi/doi/10.1073/pnas.1916387117)


Author: Neus Sabater (neus.sabater@fmi.fi)
Date: 2024-10-10
Version: 1.0

Usage:
    python convertPeatlandTifToNC.py 
    --input 
        Histosol/Histel maps in tiff
        latlonArctic grid used in the study 

 --output 
        Histosol/Histel maps in *nc at the latlonArctic
        grid used in the study

License: open-source software

Dependencies: 
    - json 
    - os
    - rasterio 
    - numpy
    - scipy
    - matplotlib
    - netCDF4
"""

import json
import os
import rasterio
from rasterio.warp import calculate_default_transform, reproject, Resampling
import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from netCDF4 import Dataset

# modify to your own path

# Define input and output files Histel peatlands
#input_tif = 'database_input/Hugelius_etal_2020_PNAS_supplement_grids_10km_Int_WAED/Histel_SOC_hg_per_sqm.tif'
#output_tif = 'database_input/Hugelius_etal_2020_PNAS_supplement_grids_10km_Int_WAED/Histel_SOC_hg_per_sqm_reprojected.tif'
#input_netcdf = 'codes/data/Computed_data/latlon_Arctic.nc'
#output_netcdf = 'codes/data/Computed_data/Histel_SOC_hg_per_sqm_reprojected.nc'

# Define input and output files Histosol peatlands
input_tif = 'database_input/Hugelius_etal_2020_PNAS_supplement_grids_10km_Int_WAED/Histosol_SOC_hg_per_sqm.tif'
output_tif = 'database_input/Hugelius_etal_2020_PNAS_supplement_grids_10km_Int_WAED/Histosol_SOC_hg_per_sqm_reprojected.tif'
input_netcdf = 'codes/data/Computed_data/latlon_Arctic.nc'
output_netcdf = 'codes/data/Computed_data/Histosol_SOC_hg_per_sqm_reprojected.nc'

# Read the input GeoTIFF file
with rasterio.open(input_tif) as src:
    data = src.read(1)
    src_crs = src.crs
    src_transform = src.transform

    # Define the target CRS (WGS84)
    dst_crs = 'EPSG:4326'

    # Calculate the transform and dimensions for the output file
    transform, width, height = calculate_default_transform(
        src_crs, dst_crs, src.width, src.height, *src.bounds)

    # Create an empty array for the reprojected data
    dst_data = np.empty((height, width), dtype=rasterio.float32)

    # Perform the reprojection
    with rasterio.open(output_tif, 'w', driver='GTiff', height=height, width=width,
                       count=1, dtype=dst_data.dtype, crs=dst_crs, transform=transform) as dst:
        reproject(
            source=rasterio.band(src, 1),
            destination=dst_data,
            src_transform=src_transform,
            src_crs=src_crs,
            dst_transform=transform,
            dst_crs=dst_crs,
            resampling=Resampling.nearest)
        dst.write(dst_data, 1)

# Read the latitude and longitude vectors from the netCDF file
with Dataset(input_netcdf, 'r') as nc:
    lat_vector = nc.variables['latitude'][:]
    lon_vector = nc.variables['longitude'][:]

# Extract the coordinates and values from the reprojected data
with rasterio.open(output_tif) as src:
    reprojected_data = src.read(1)
    transform = src.transform
    
    # Generate the meshgrid of pixel coordinates
    x_coords, y_coords = np.meshgrid(np.arange(reprojected_data.shape[1]), np.arange(reprojected_data.shape[0]))
    
    # Convert pixel coordinates to geographic coordinates (lon, lat)
    lon, lat = rasterio.transform.xy(transform, y_coords.flatten(), x_coords.flatten(), offset='center')
    
    # Convert lon and lat lists to numpy arrays
    lon = np.array(lon)
    lat = np.array(lat)

    # Reshape lon and lat to match the shape of reprojected_data
    lon = lon.reshape(reprojected_data.shape)
    lat = lat.reshape(reprojected_data.shape)


# Flatten the arrays for interpolation
lon_flat = lon.flatten()
lat_flat = lat.flatten()
values_flat = reprojected_data.flatten()

# Generate the lat/lon grid from the vectors
lon_grid, lat_grid = np.meshgrid(lon_vector, lat_vector)

# Interpolate the reprojected data to the new grid
new_data = griddata((lon_flat, lat_flat), values_flat, (lon_grid, lat_grid), method='linear')

# Plot the regridded data
plt.figure()
plt.imshow(new_data, extent=(lon_vector.min(), lon_vector.max(), lat_vector.min(), lat_vector.max()), origin='lower')
plt.colorbar(label='Value')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.title('Reprojected and Regridded Data')
plt.show()

# Save the regridded data to a new netCDF file
with Dataset(output_netcdf, 'w', format='NETCDF4') as nc_out:
    # Create dimensions
    nc_out.createDimension('lat', len(lat_vector))
    nc_out.createDimension('lon', len(lon_vector))
    
    # Create variables
    latitudes = nc_out.createVariable('lat', 'f4', ('lat',))
    longitudes = nc_out.createVariable('lon', 'f4', ('lon',))
    regridded_data = nc_out.createVariable('regridded_data', 'f4', ('lat', 'lon'))
    
    # Assign data to variables
    latitudes[:] = lat_vector
    longitudes[:] = lon_vector
    regridded_data[:, :] = new_data
    
    # Add attributes
    latitudes.units = 'degrees_north'
    longitudes.units = 'degrees_east'
    regridded_data.units = 'hg / m-2'
    #regridded_data.long_name = 'Estimated storage of organic carbon in peat of Histel peatlands (unit is hg / m-2)'
    regridded_data.long_name = 'Estimated storage of organic carbon in peat of Histosol peatlands (unit is hg / m-2)'


print(f'Regridded data saved to {output_netcdf}')
