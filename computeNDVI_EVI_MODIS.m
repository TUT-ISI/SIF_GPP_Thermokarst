% ###################################################################
% ####### computeNDVI_EVI_MODIS ####
% This script reads the EVI and NDVI products from MODIS Terra and Aqua
% satellites and regrid them in the ArcticSIF 0.1 degree resolution and in
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

path_in  = config.input_dir.NDVI; % Here there are two folders from AQUA and TERRA
path_out = config.input_dir.computed; % Saving the results in the computed data folder for a later use 

saving_flag   = 1;
plotting_flag = 0;
skipping_initial = 0; % Skip the reading from MODIS data and directly execute the merging between AQUAN and TERRA extracted files


filesAndFolders = dir(path_in);
% Extract the names and whether each item is a directory
names = {filesAndFolders.name};
isDir = [filesAndFolders.isdir];

folders = names(isDir);
folders = folders(~ismember(folders, {'.', '..'}));
num_dir  = length(folders);

% load the info regarding the ArcticSIF geometry coverage (60-90)N
%ncdisp([config.input_dir.GEO,'latlon_Arctic.nc'])
latitude  = ncread([config.input_dir.GEO,'latlon_Arctic.nc'],'latitude');
longitude = ncread([config.input_dir.GEO,'latlon_Arctic.nc'],'longitude');

lat_MODIS = [89.975:-0.05:-89.975];
lon_MODIS = [-179.975:0.05:179.975];

% Preparing data to reproject to the Arctic region used in ArcticSIF 60 N - 90 N
[LAT,LON] = meshgrid(latitude,longitude);
[LAT_MODIS, LON_MODIS] = meshgrid(lat_MODIS,lon_MODIS);
LAT_MODIS         = LAT_MODIS';
LON_MODIS         = LON_MODIS';

% Flatten the MODIS grids into 1D vectors
LAT_MODIS_flat = LAT_MODIS(:);
LON_MODIS_flat = LON_MODIS(:);

% Precompute the nearest neighbor indices to reproject faster
% F = scatteredInterpolant(LAT_MODIS_flat, LON_MODIS_flat, (1:numel(LAT_MODIS_flat))', 'nearest', 'nearest');
% nearest_indices = F(LAT, LON);
% nearest_indices = single(nearest_indices);
load([config.input_dir.IGBP,filesep,'nearest_indices_MODIS_ArcticSIF.mat'])


% Land cover classification  % '16 - Water Bodies
% ---------------------------
% path_igbp           = config.input_dir.IGBP;
% igbp_RE             = ncread([path_igbp,filesep,'IGBP_05d.nc'],'IGBP');
% lat_igbp            = ncread([path_igbp,filesep,'IGBP_05d.nc'],'lat_igbp');
% lon_igbp            = ncread([path_igbp,filesep,'IGBP_05d.nc'],'lon_igbp');
% class_val           = unique(igbp_RE);
% num_class           = length(class_val);
% [LAT_igbp,LON_igbp] = meshgrid(lat_igbp,lon_igbp);
% LAT_igbp_flat       = double(LAT_igbp(:));
% LON_igbp_flat       = double(LON_igbp(:));
% igbp_RE_flat        = double(igbp_RE(:));
% F                   = scatteredInterpolant(LAT_igbp_flat, LON_igbp_flat, double(igbp_RE_flat), 'nearest', 'none'); % Create a scatteredInterpolant object for nearest neighbor interpolation
% % Interpolate the data on the new grid
% igbp_Arctic         = F(LAT, LON);
% clear lat_igbp lon_igbp LAT_igbp LON_igbp igbp_RE
% land_mask = ones(size(igbp_Arctic));
% land_mask(igbp_Arctic==16) = NaN;
% save([config.input_dir.IGBP,filesep,'land_mask.mat'],'land_mask')
load([config.input_dir.IGBP,filesep,'land_mask.mat'])


% Hard coded info from the NDVI and EVI files 8read)
VI_info.valid_range = [0 65534];
EVI_info.scalefactor = 10000;
EVI_info.fillValue = -3000;
EVI_info.offset = 0;
EVI_info.valid_range = [-2000,10000];

NDVI_info.scalefactor = 10000;
NDVI_info.fillValue = -3000;
NDVI_info.offset = 0;
NDVI_info.valid_range = [-2000,10000];


for i=1:num_dir
    files = dir([path_in,filesep,folders{i}]);
    files = files(3:end); % First two are {'.', '..'}

    % pre-allocate memory 
    % JUL_DAY           = NaN(1,length(files));
    % MONTH_            = NaN(1,length(files));
    % DAY_              = NaN(1,length(files));
    % YEAR_             = NaN(1,length(files));
    % NDVI_Arctic       = NaN(length(files),length(longitude),length(latitude));
    % EVI_Arctic        = NaN(length(files),length(longitude),length(latitude));
    % NDVI_STD_Arctic   = NaN(length(files),length(longitude),length(latitude));
    % EVI_STD_Arctic    = NaN(length(files),length(longitude),length(latitude));
    % VI_QUALITY_Arctic = NaN(length(files),length(longitude),length(latitude));
    j = 0;
    for k = 1:length(files)
        disp(['folder',num2str(i),'- file', num2str(k)])

        filename = [path_in,filesep,folders{i},filesep,files(k).name];
        % Extract the year and Julian day
        aux       = files(k).name;
        year      = str2double(aux(10:13));
        julianDay = str2double(aux(14:16));

        % Convert Julian day to calendar date
        dateNum      = datenum(year, 1, 1) + julianDay - 1;
        calendarDate = datestr(dateNum, 'yyyy-mm-dd');
        month_       = str2num(calendarDate(6:7));

        if (month_==6 || month_ == 7) % Store only summer months
            j = j + 1;             

            info                = hdfinfo([path_in,filesep,folders{i},filesep,files(k).name]);

            NDVI                = hdfread(filename,'/MODIS_Grid_16Day_VI_CMG/Data Fields/CMG 0.05 Deg 16 days NDVI', 'Index', {[1  1],[1  1],[3600  7200]});
            pixel_reliability   = hdfread(filename,'/MODIS_Grid_16Day_VI_CMG/Data Fields/CMG 0.05 Deg 16 days pixel reliability', 'Index', {[1  1],[1  1],[3600  7200]});

            EVI                 = hdfread(filename, '/MODIS_Grid_16Day_VI_CMG/Data Fields/CMG 0.05 Deg 16 days EVI', 'Index', {[1  1],[1  1],[3600  7200]});
            VI_Quality          = hdfread(filename, '/MODIS_Grid_16Day_VI_CMG/Data Fields/CMG 0.05 Deg 16 days VI Quality', 'Index', {[1  1],[1  1],[3600  7200]});

            NDVI_std_dev        = hdfread(filename, '/MODIS_Grid_16Day_VI_CMG/Data Fields/CMG 0.05 Deg 16 days NDVI std dev', 'Index', {[1  1],[1  1],[3600  7200]});
            EVI_std_dev         = hdfread(filename, '/MODIS_Grid_16Day_VI_CMG/Data Fields/CMG 0.05 Deg 16 days EVI std dev', 'Index', {[1  1],[1  1],[3600  7200]});


            ndvi                = double(NDVI)./NDVI_info.scalefactor;
            evi                 = double(EVI)./ EVI_info.scalefactor;
            ndvi_std            = double(NDVI_std_dev)./NDVI_info.scalefactor;
            evi_std             = double(EVI_std_dev)./ EVI_info.scalefactor;



            % Reproject variables A, B, C using the precomputed indices
            ndvi_Arctic         = reproject_variable(ndvi, nearest_indices);
            evi_Arctic          = reproject_variable(evi, nearest_indices);
            ndvi_std_Arctic     = reproject_variable(ndvi_std, nearest_indices);
            evi_std_Arctic      = reproject_variable(evi_std, nearest_indices);
            VI_Quality_Arctic   = reproject_variable(VI_Quality, nearest_indices);

            % Reshape the reprojected variables to match the new grid shape
            ndvi_Arctic         = reshape(ndvi_Arctic, size(LAT));
            evi_Arctic          = reshape(evi_Arctic, size(LAT));
            ndvi_std_Arctic     = reshape(ndvi_std_Arctic, size(LAT));
            evi_std_Arctic      = reshape(evi_std_Arctic, size(LAT));
            VI_Quality_Arctic   = reshape(VI_Quality_Arctic, size(LAT));

            % Put NaNs (exclude) non-land surface uisng the IGBP
            % classification
            ndvi_Arctic         = ndvi_Arctic.*land_mask;
            evi_Arctic          = evi_Arctic.*land_mask;
            ndvi_std_Arctic     = ndvi_std_Arctic.*land_mask;
            evi_std_Arctic      = evi_std_Arctic.*land_mask;



            % Store the data in ram memory
            JUL_DAY(j)               = julianDay;
            CALENDAR{j}              = calendarDate;
            YEAR_(j)                 = str2num(calendarDate(1:4));
            MONTH_(j)                = str2num(calendarDate(6:7));
            DAY_(j)                  = str2num(calendarDate(9:10));
            NDVI_Arctic(j,:,:)       = ndvi_Arctic;
            EVI_Arctic(j,:,:)        = evi_Arctic;
            NDVI_STD_Arctic(j,:,:)   = ndvi_std_Arctic;
            EVI_STD_Arctic(j,:,:)    = evi_std_Arctic;
            VI_QUALITY_Arctic(j,:,:) = VI_Quality_Arctic;

        end

    end
    num_steps = j; % Total number of MODIS data in June-July

    % Store the data in a separate NetCDF file for AQUA and TERRA

    % NetCDF file
    filename = [path_out,folders{i},'.nc'];

    % Create a new NetCDF file format 4 to allocate more memory
    ncid = netcdf.create(filename, 'NETCDF4');

    % Define dimensions
    dimid_lon = netcdf.defDim(ncid, 'lon', length(longitude));
    dimid_lat = netcdf.defDim(ncid, 'lat', length(latitude));
    dimid_t   = netcdf.defDim(ncid, 'time_steps', num_steps);

    % Define variables
    varid_NDVI              = netcdf.defVar(ncid, 'NDVI_Arctic', 'NC_FLOAT', [dimid_t,dimid_lon,dimid_lat]);
    varid_EVI               = netcdf.defVar(ncid, 'EVI_Arctic', 'NC_FLOAT', [dimid_t,dimid_lon,dimid_lat]);
    varid_NDVI_STD_Arctic   = netcdf.defVar(ncid, 'NDVI_STD_Arctic', 'NC_FLOAT', [dimid_t,dimid_lon,dimid_lat]);
    varid_EVI_STD_Arctic    = netcdf.defVar(ncid, 'EVI_STD_Arctic', 'NC_FLOAT', [dimid_t,dimid_lon,dimid_lat]);
    varid_VI_QUALITY_Arctic = netcdf.defVar(ncid, 'VI_QUALITY_Arctic', 'NC_FLOAT', [dimid_t,dimid_lon,dimid_lat]);

    varid_latitude          = netcdf.defVar(ncid, 'latitude', 'NC_FLOAT', dimid_lat);
    varid_longitude         = netcdf.defVar(ncid, 'longitude', 'NC_FLOAT', dimid_lon);
    varid_JUL_DAY           = netcdf.defVar(ncid, 'JUL_DAY', 'NC_FLOAT', dimid_t);
    varid_MONTH             = netcdf.defVar(ncid, 'MONTH', 'NC_FLOAT', dimid_t);
    varid_DAY               = netcdf.defVar(ncid, 'DAY', 'NC_FLOAT', dimid_t);
    varid_YEAR              = netcdf.defVar(ncid, 'YEAR', 'NC_FLOAT', dimid_t);



    % Add attributes for descriptions
    netcdf.putAtt(ncid, varid_NDVI, 'description', 'NDVI values for the ArcticSIF region MYD13C1 v006 MODIS product');
    netcdf.putAtt(ncid, varid_EVI, 'description', 'EVI values for the ArcticSIF region MYD13C1 v006 MODIS product');
    netcdf.putAtt(ncid, varid_NDVI_STD_Arctic, 'description', 'NDVI std values for the ArcticSIF region MYD13C1 v006 MODIS product');
    netcdf.putAtt(ncid, varid_EVI_STD_Arctic, 'description', 'EVI std values for the ArcticSIF region MYD13C1 v006 MODIS product');
    netcdf.putAtt(ncid, varid_VI_QUALITY_Arctic, 'description',['Quality index values for the ArcticSIF region MYD13C1 v006 MODIS ',...
                                                ' Rank Keys: [-1]: Fill/No Data -Not Processed [0]: Good data - Use with confidence',...
                                                '[1]: Marginal data - Useful, but look at other QA information',...
                                                '[2]: Snow/Ice - Target covered with snow/ice',...
                                                '[3]: Cloudy - Target not visible, covered with cloud',...
                                                '[4]: Estimated - From MODIS historic time series']);

    netcdf.putAtt(ncid, varid_JUL_DAY, 'description', 'Julian day of the year for the ArcticSIF region MYD13C1 v006 MODIS product');
    netcdf.putAtt(ncid, varid_YEAR, 'description', 'Year for the ArcticSIF region MYD13C1 v006 MODIS product');
    netcdf.putAtt(ncid, varid_MONTH, 'description', 'Month for the ArcticSIF region MYD13C1 v006 MODIS product');
    netcdf.putAtt(ncid, varid_DAY, 'description', 'Day for the ArcticSIF region MYD13C1 v006 MODIS product');

    
    % End definitions and enter data mode
    netcdf.endDef(ncid);


    % Store the data into the NetCDF variables
    netcdf.putVar(ncid, varid_NDVI, single(NDVI_Arctic));
    netcdf.putVar(ncid, varid_EVI, single(EVI_Arctic));
    netcdf.putVar(ncid, varid_NDVI_STD_Arctic, single(NDVI_STD_Arctic));
    netcdf.putVar(ncid, varid_EVI_STD_Arctic,single( EVI_STD_Arctic));
    netcdf.putVar(ncid, varid_VI_QUALITY_Arctic, single(VI_QUALITY_Arctic));
    netcdf.putVar(ncid, varid_latitude, single(latitude));
    netcdf.putVar(ncid, varid_longitude, single(longitude));
    netcdf.putVar(ncid, varid_JUL_DAY,single(JUL_DAY));
    netcdf.putVar(ncid, varid_YEAR, single(YEAR_));
    netcdf.putVar(ncid, varid_MONTH, single(MONTH_));
    netcdf.putVar(ncid, varid_DAY, single(DAY_));



    % Close the NetCDF file
    netcdf.close(ncid);
    disp(['Saved NetCDF file'])

end


% Merging TERRA and AQUA datasets done in
% computeNDVI_EVI_MODISAQUATERRA_Merging.m



% Store the full data in a NetCDF:
% -NDVI_JJ [lat,lon,years]
% -NDVI_MA [lat,lon,years]
% -NDVI_AS [lat,lon,years]
% The same for EVI and the same for the std deviation





