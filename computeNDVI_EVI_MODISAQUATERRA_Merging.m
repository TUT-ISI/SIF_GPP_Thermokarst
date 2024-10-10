% ###################################################################
% ####### computeNDVI_EVI_MODISAQUATERRA_Merging ####
% This script reads the EVI and NDVI NetCDF products from MODIS Terra and Aqua
% satellites and regrid then in the ArcticSIF 0.1 degree resolution and in
% an 8-day interval. 

% Author: Neus Sabater
% e-mail: neus.sabater@fmi.fi 
% version 0: (28/05/2024) 
% ------------------------------------------------------
% ------------------------------------------------------

% Main script
addpath(genpath([pwd,filesep,'functions_environment']))
addpath(genpath([pwd,filesep,'functions_matlab']))

% Automatically detect environment
env = detect_environment();

% Load configuration
config = load_conf(env);

path_in  = config.input_dir.computed; % Saving the previous AQUA and TERRA results in the computed data folder
path_out = config.input_dir.computed; % Saving the results in the computed data folder for a later use 


% Read the NetCDF containing the datetimes from AQUA in
% June/July
aqua_file   = [path_in,'NDVI_EVI_AQUA.nc'];
jul_aqua    = ncread(aqua_file,'JUL_DAY'); 
month_aqua  = ncread(aqua_file,'MONTH');
day_aqua    = ncread(aqua_file,'DAY');
year_aqua   = ncread(aqua_file,'YEAR');

latitude    = ncread(aqua_file,'latitude');% same for Aqua and Terra
longitude   = ncread(aqua_file,'longitude');% same for Aqua and Terra
info        = ncinfo(aqua_file);
dimensions  = info.Dimensions;
for i = 1:length(dimensions)
    eval([dimensions(i).Name,'= dimensions(i).Length;']) % same for Aqua and Terra
end

% Read the NetCDF containing the datetimes from TERRA in
% June/July
terra_file   = [path_in,'NDVI_EVI_TERRA.nc'];
jul_terra    = ncread(terra_file,'JUL_DAY'); 
month_terra  = ncread(terra_file,'MONTH');
day_terra    = ncread(terra_file,'DAY');
year_terra   = ncread(terra_file,'YEAR');


%Group per year: Jun/Jul AQUA and TERRA year 20xx
year_list = 2003:2018;
for i=1:length(year_list)
    year_list(i)

    [pos_terra,~]   = find(year_terra==year_list(i));
    [pos_aqua,~]    = find(year_aqua==year_list(i));

    for j = 1:length(pos_aqua)
        ndvi_aqua(j,:,:)       = ncread(aqua_file,'NDVI_Arctic',[pos_aqua(j),1,1],[1,lon,lat]);
        evi_aqua(j,:,:)        = ncread(aqua_file,'EVI_Arctic',[pos_aqua(j),1,1],[1,lon,lat]);
    end
    for j = 1:length(pos_terra)
        ndvi_terra(j,:,:)       = ncread(terra_file,'NDVI_Arctic',[pos_terra(j),1,1],[1,lon,lat]);
        evi_terra(j,:,:)        = ncread(terra_file,'EVI_Arctic',[pos_terra(j),1,1],[1,lon,lat]);
    end
    
    ndvi_aux        = [ndvi_aqua;ndvi_terra];
    ndvi_all        = mean(ndvi_aux,1,'omitnan');
    ndvi_std        = std(ndvi_aux,1,'omitnan');

    NDVI(i,:,:)     = ndvi_all;
    NDVI_std(i,:,:) = ndvi_std;

    evi_aux        = [evi_aqua;evi_terra];
    evi_all        = mean(evi_aux,1,'omitnan');
    evi_std        = std(evi_aux,1,'omitnan');

    EVI(i,:,:)     = evi_all;
    EVI_std(i,:,:) = evi_std;

end


%store in a NetCDF file
% NetCDF file
filename = [path_out,'NDVI_EVI_yearly.nc'];

% Create a new NetCDF file format 4 to allocate more memory
ncid = netcdf.create(filename, 'NETCDF4');

% Define dimensions
dimid_lon = netcdf.defDim(ncid, 'lon', lon);
dimid_lat = netcdf.defDim(ncid, 'lat', lat);
num_steps = length(year_list);
dimid_t   = netcdf.defDim(ncid, 'time_steps', num_steps);

% Define variables
varid_NDVI              = netcdf.defVar(ncid, 'NDVI_Arctic', 'NC_FLOAT', [dimid_t,dimid_lon,dimid_lat]);
varid_EVI               = netcdf.defVar(ncid, 'EVI_Arctic', 'NC_FLOAT', [dimid_t,dimid_lon,dimid_lat]);
varid_NDVI_STD_Arctic   = netcdf.defVar(ncid, 'NDVI_STD_Arctic', 'NC_FLOAT', [dimid_t,dimid_lon,dimid_lat]);
varid_EVI_STD_Arctic    = netcdf.defVar(ncid, 'EVI_STD_Arctic', 'NC_FLOAT', [dimid_t,dimid_lon,dimid_lat]);


varid_latitude          = netcdf.defVar(ncid, 'latitude', 'NC_FLOAT', dimid_lat);
varid_longitude         = netcdf.defVar(ncid, 'longitude', 'NC_FLOAT', dimid_lon);
varid_YEAR              = netcdf.defVar(ncid, 'YEAR', 'NC_FLOAT', dimid_t);


% Add attributes for descriptions
netcdf.putAtt(ncid, varid_NDVI, 'description', 'Average NDVI values for the ArcticSIF region MYD13C1 v006 MODIS product TERRA and AQUA merged');
netcdf.putAtt(ncid, varid_EVI, 'description', 'Average EVI  values for the ArcticSIF region MYD13C1 v006 MODIS product TERRA and AQUA merged');
netcdf.putAtt(ncid, varid_NDVI_STD_Arctic, 'description', 'Standard deviation of the NDVI values for the ArcticSIF region MYD13C1 v006 MODIS product TERRA and AQUA merged');
netcdf.putAtt(ncid, varid_EVI_STD_Arctic, 'description', 'Standard deviation of the EVI values for the ArcticSIF region MYD13C1 v006 MODIS product TERRA and AQUA merged');
netcdf.putAtt(ncid, varid_YEAR, 'description', 'Year for the ArcticSIF region MYD13C1 v006 MODIS product');


% End definitions and enter data mode
netcdf.endDef(ncid);


% Store the data into the NetCDF variables
netcdf.putVar(ncid, varid_NDVI, single(NDVI));
netcdf.putVar(ncid, varid_EVI, single(EVI));
netcdf.putVar(ncid, varid_NDVI_STD_Arctic, single(NDVI_std));
netcdf.putVar(ncid, varid_EVI_STD_Arctic,single(EVI_std));
netcdf.putVar(ncid, varid_latitude, single(latitude));
netcdf.putVar(ncid, varid_longitude, single(longitude));
netcdf.putVar(ncid, varid_YEAR, single(year_list));

% Close the NetCDF file
netcdf.close(ncid);
disp(['Saved NetCDF file'])





