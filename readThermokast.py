"""
Script: Extract and Project Data from Shapefile of the Thermokast database
Author: Neus Sabater
Version: 1.0
Date: March 29, 2024

Description:
This script opens a Shapefile dataset containing fields: TKWP, TKThLP, TKHP, and TSOC_kgC,
along with their corresponding coordinates. It extracts these fields, projects them onto a
latitude and longitude grid, and saves the data in a NetCDF file named 'latlon_Arctic.nc'.
Paper available in 
https://www.nature.com/articles/ncomms13043
Data available in 
https://daac.ornl.gov/cgi-bin/dsviewer.pl?ds_id=1332

"""

import geopandas as gpd
import xarray as xr
from pyproj import Proj, transform
from scipy.interpolate import griddata
import os
import numpy as np
from netCDF4 import Dataset
from sklearn.preprocessing import OrdinalEncoder
from scipy.spatial import cKDTree

# Function to read shapefile and extract data
def read_shapefile(shapefile):
    # Read the shapefile
    gdf = gpd.read_file(shapefile)

    # Extract desired fields
    tkwp = gdf['TKWP']
    tkthlp = gdf['TKThLP']
    tkhp = gdf['TKHP']
    tsoc_kgc = gdf['TSOC_kgC']

    # Get coordinates (centroid of each polygon)
    gdf['lon'] = gdf.centroid.x
    gdf['lat'] = gdf.centroid.y
    lon, lat = gdf['lon'].values, gdf['lat'].values

    return tkwp, tkthlp, tkhp, tsoc_kgc, lon, lat

# Function to convert coordinates to lat/lon
def convert_to_latlon(x, y, proj4_str):
    input_proj = Proj(proj4_str)
    output_proj = Proj(init='epsg:4326')  # WGS84, standard lat/lon

    lon, lat = transform(input_proj, output_proj, x, y)
    return lat, lon

# Function to create netCDF file
def create_netcdf(tkwp, tkthlp, tkhp, tsoc_kgc, lat, lon, output_file):
    # Create xarray dataset
    ds = xr.Dataset({
        'TKWP': (('lat', 'lon'), tkwp),
        'TKThLP': (('lat', 'lon'), tkthlp),
        'TKHP': (('lat', 'lon'), tkhp),
        'TSOC_kgC': (('lat', 'lon'), tsoc_kgc),
    },
    coords={'lat': lat, 'lon': lon})

    # Save to netCDF file
    ds.to_netcdf(output_file)
    print("NetCDF file saved successfully:", output_file)

# Main function
def main():
    # Path to shapefile and projection file
    path_in   = '/Users/neus/Dropbox/IPL/02_PROYECTOS_PROPIOS/SIF_GPP/database_input/PERMAFROST/Thermokarst_Circumpolar_Map_1332/data/Circumpolar_Thermokarst_Landscapes/'
    shapefile = os.path.join(path_in, 'Circumpolar_Thermokarst_Landscapes.shp')
    proj_file = os.path.join(path_in, 'Circumpolar_Thermokarst_Landscapes.prj')
    path_out  = '/Users/neus/Dropbox/IPL/02_PROYECTOS_PROPIOS/SIF_GPP/codes/data/Computed_data/'
    geofile   = os.path.join(path_out,'latlon_Arctic.nc')
    
    # Read latitude and longitude fields from the NetCDF file
    try:
    # Open the NetCDF file
        ds = xr.open_dataset(geofile)

        # Read latitude and longitude arrays
        latitude = ds['latitude']
        longitude = ds['longitude']

        # Print the first few values as an example
        print("Latitude values:")
        print(latitude.values[:10])  # Print first 10 values
        print("Longitude values:")
        print(longitude.values[:10])  # Print first 10 values

        # Close the NetCDF file
        ds.close()

    except Exception as e:
        print(f"Error reading NetCDF file: {e}")

    # Read shapefile and extract data
    tkwp, tkthlp, tkhp, tsoc_kgc, x, y = read_shapefile(shapefile)

    # Extract projection information
    with open(proj_file, 'r') as f:
        proj4_str = f.read()

    # Convert coordinates to lat/lon
    lat, lon = convert_to_latlon(x, y, proj4_str)

    # Create meshgrid from latitude and longitude
    lon_mesh, lat_mesh = np.meshgrid(longitude, latitude)

    # # Categorical data
    # categories = ['None','Low','Moderate','High','Very High']
    # numeric_values = [0,1,2,3,4]


    # #  Handling missing data and replace empty strings with a placeholder, such as 'Missing'
    # tkwp = ['Missing' if x == '' else x for x in tkwp]
    # tkthlp= ['Missing' if x == '' else x for x in tkthlp]
    # tkhp = ['Missing' if x == '' else x for x in tkhp]

  
    cat2num = {'None': 0, 'Low': 1, 'Moderate': 2, 'High': 3, 'Very High': 4, None: 5}
    tkwp_encoded = np.array([cat2num[num] for num in tkwp])
    tkthlp_encoded =  np.array([cat2num[num] for num in tkthlp])
    tkhp_encoded =  np.array([cat2num[num] for num in tkhp])

  
    # Check for NaN values in tsoc_kgc, lon, and lat and remove them
    nan_mask = np.isnan(tsoc_kgc) | np.isnan(lon) | np.isnan(lat)
    valid_indices = ~nan_mask
    lon_valid = lon[valid_indices]
    lat_valid = lat[valid_indices]
    tsoc_kgc_valid = tsoc_kgc[valid_indices]

    # Interpolate data to latitude and longitude grid without NaNs
    tsoc_kgc_regrid = griddata((lon_valid, lat_valid), tsoc_kgc_valid, (lon_mesh, lat_mesh), method='linear',fill_value=np.nan)

    # Calculate distances between interpolated points and nearest points in original data
    tree = cKDTree(np.column_stack((lon_valid, lat_valid)))
    distances, indices = tree.query(np.column_stack((lon_mesh.flatten(), lat_mesh.flatten())))

    # Convert distances from flat array to 2D grid
    distances_grid = distances.reshape(lon_mesh.shape)

    # Set threshold distance in latitude and longitude
    threshold_lat = 0.3  # Example threshold distance in latitude (adjust as needed)
    threshold_lon = 0.3  # Example threshold distance in longitude (adjust as needed)
    # Set points in interpolated grid to NaN where distance exceeds threshold
    tsoc_kgc_regrid[distances_grid > threshold_lon] = np.nan
    tsoc_kgc_regrid[distances_grid > threshold_lat] = np.nan

    # Check for NaN values in tkwp_encoded, lon, and lat and remove them
    nan_mask = np.isnan(tkwp_encoded) | np.isnan(lon) | np.isnan(lat)
    valid_indices = ~nan_mask
    lon_valid = lon[valid_indices]
    lat_valid = lat[valid_indices]
    tkwp_encoded_valid = tkwp_encoded[valid_indices]
    tkwp_regrid = griddata((lon_valid, lat_valid), tkwp_encoded_valid, (lon_mesh, lat_mesh), method='nearest',fill_value=np.nan)
    tkwp_regrid = tkwp_regrid.astype(float)
    tkwp_regrid[distances_grid > threshold_lon] = np.nan
    tkwp_regrid[distances_grid > threshold_lat] = np.nan

    # Check for NaN values in tkthlp_encoded, lon, and lat and remove them
    nan_mask = np.isnan(tkthlp_encoded) | np.isnan(lon) | np.isnan(lat)
    valid_indices = ~nan_mask
    lon_valid = lon[valid_indices]
    lat_valid = lat[valid_indices]
    tkthlp_encoded_valid = tkthlp_encoded[valid_indices]
    tkthlp_regrid = griddata((lon_valid, lat_valid), tkthlp_encoded_valid, (lon_mesh, lat_mesh), method='nearest',fill_value=np.nan)
    tkthlp_regrid = tkthlp_regrid.astype(float)
    tkthlp_regrid[distances_grid > threshold_lon] = np.nan
    tkthlp_regrid[distances_grid > threshold_lat] = np.nan

    # Check for NaN values in tkhp_encoded, lon, and lat and remove them
    nan_mask = np.isnan(tkhp_encoded) | np.isnan(lon) | np.isnan(lat)
    valid_indices = ~nan_mask
    lon_valid = lon[valid_indices]
    lat_valid = lat[valid_indices]
    tkhp_encoded_valid = tkhp_encoded[valid_indices]
    tkhp_regrid = griddata((lon_valid, lat_valid), tkhp_encoded_valid, (lon_mesh, lat_mesh), method='nearest',fill_value=np.nan)
    tkhp_regrid = tkhp_regrid.astype(float)
    tkhp_regrid[distances_grid > threshold_lon] = np.nan
    tkhp_regrid[distances_grid > threshold_lat] = np.nan

    # Combine the three thermokast datasets in only one. 
    # Applying the criterium that the highest values of all of three.
    # create aux file with the 5 values (missing data) replaced to -1
    tkwp_aux = tkwp_regrid
    tkwp_aux[tkwp_aux==5] = -1
    tkthlp_aux = tkthlp_regrid
    tkthlp_aux[tkthlp_aux==5] = -1
    tkhp_aux = tkhp_regrid
    tkhp_aux[tkhp_aux==5] = -1
    
    tk_all_regrid =  np.full_like(tkwp_regrid, np.nan)
    # Iterate through each element of the matrices and take the maximum
    for i in range(tkwp_aux.shape[0]):
        for j in range(tkwp_aux.shape[1]):
            tk_all_regrid[i, j] = max(tkwp_aux[i, j], tkthlp_aux[i, j], tkhp_aux[i, j])


    output_file = os.path.join(path_out,'Thermokast_info_all.nc')
    # Create netCDF file
    # create_netcdf(tkwp_regrid, tkthlp_regrid, tkhp_regrid, tsoc_kgc_regrid, latitude, longitude, output_file)
    
    ncout = Dataset(output_file, 'w', format='NETCDF4')

    # Add some attributes
    ncout.History = 'File generated from the Thermokast information https://daac.ornl.gov/cgi-bin/dsviewer.pl?ds_id=1332'
    ncout.Description = 'Thermokast on wetland, hill and lake sides. Also total carbon soil Carbon info.'

    # create dimensions
    ncout.createDimension('latitude', len(latitude))
    ncout.createDimension('longitude', len(longitude))

    # save coordinates
    ncout_latitude = ncout.createVariable('latitude', 'f4', ('latitude',))
    ncout_latitude[:] = latitude
    ncout_longitude = ncout.createVariable('longitude', 'f4', ('longitude',))
    ncout_longitude[:] = longitude

    # save data
    ncout_dist = ncout.createVariable('tkwp', 'f4', ('latitude', 'longitude'))
    ncout_dist.long_name = 'Wetland thermokarst terrain coverage'
    ncout_dist.units = 'None, Low, Moderate, High, and Very High expressed as numbers from 0 to 4 '
    ncout_dist[:] = tkwp_regrid

    ncout_dist = ncout.createVariable('tkthlp', 'f4', ('latitude', 'longitude'))
    ncout_dist.long_name = 'Lake thermokarst terrain coverage'
    ncout_dist.units = 'None, Low, Moderate, High, and Very High expressed as numbers from 0 to 4 '
    ncout_dist[:] = tkthlp_regrid

    ncout_dist = ncout.createVariable('tkhp', 'f4', ('latitude', 'longitude'))
    ncout_dist.long_name = 'Hillslope thermokarst terrain coverage'
    ncout_dist.units = 'None, Low, Moderate, High, and Very High expressed as numbers from 0 to 4 '
    ncout_dist[:] = tkhp_regrid

    ncout_dist = ncout.createVariable('tk_all', 'f4', ('latitude', 'longitude'))
    ncout_dist.long_name = 'Max thermokarst terrain coverage including wetland, lake and hillslope'
    ncout_dist.units = 'None, Low, Moderate, High, and Very High expressed as numbers from 0 to 4 '
    ncout_dist[:] = tk_all_regrid

    ncout_dist = ncout.createVariable('TSOC_kgC', 'f4', ('latitude', 'longitude'))
    ncout_dist.long_name = 'Total Soil Organic Carbon in Kg C'
    ncout_dist.units = ' Kg C'
    ncout_dist[:] = tsoc_kgc_regrid


    ncout.close()
    print('Done!')

if __name__ == "__main__":
    main()
