% ###################################################################
% ####### READ THE ARCTIC DATA CENTER MAP FOR CLASSIFICATION  ####
% This script read the *.shp file and *prj files for the land and
% lake arctic database published here
% https://arcticdata.io/catalog/view/doi:10.18739/A2C824F9X
% 
% Author: Neus Sabater
% e-mail: neus.sabater@fmi.fi 
% version 0: (19/02/2024) 
% ------------------------------------------------------
% ------------------------------------------------------

addpath(genpath([pwd,filesep,'functions_environment']))
% Automatically detect environment
env = detect_environment();
% Load configuration
config = load_conf(env);

% setting for saving and ploting 
path_in     = config.input_dir.computed;  % SIF GPP data (needed for the lat/lon)
path_fig    = config.output_dir.general;
save_flag   = 0;


% Specify the path to your Shapefile and its associated projection file
shpFile = [config.input_dir.general,'/arctic_data_center/resource_map_doi_10_18739_A2C824F9X/data/BAWLD_V1___Shapefile/BAWLD_V1.shp'];
prjFile = [config.input_dir.general,'/arctic_data_center/resource_map_doi_10_18739_A2C824F9X/data/BAWLD_V1___Shapefile/BAWLD_V1.prj'];

% Read the Shapefile using shaperead
S = shaperead(shpFile);


%% mapping in polar plots
% ---------------------------
addpath(genpath(config.matlab_paths.mappingToolbox))

% Open a new figure
num = length(S);

for i= 1:num
    lat (i)=S(i).Lat1;
    lon (i)=S(i).Long;
    GLA_(i)=S(i).GLA; % Percent of galciers and permanent snow within grid cell (%).  
    ROC_(i)=S(i).ROC; % Percent of rockland within grid cell (%).
    TUN_(i)=S(i).TUN; % Percent of tundra within grid cell (%).
    BOR_(i)=S(i).BOR; % Percent of boreal forest within grid cell (%).
    WET_(i)=S(i).WET; % Percent of wetland within grid cell (%) Includes, permafrost bogs, wet tundra, marsh, bog, and fen.
    LAK_(i)=S(i).LAK; % Percent of lakes within grid cell (%).
    RIV_(i)=S(i).RIV; % Percent of rivers within grid cell (%).
    CLASS_(i)=S(i).WETSCAPE; % K-means classification. 
end




%Interpolate to the lat/lon of the SIF-GPP datasets
load([path_in,'SIF_GPP_Z_score_TrendArticSIFGPPFluxSat_JJ_2003_2018.mat'],'latitude_range','longitude') 
[LAT,LON]  = meshgrid(latitude_range,longitude); % coordinates of the SIF and GPP data 
BOR_int    = griddata(lat,lon,BOR_,double(LAT),double(LON),'linear');
GLA_int    = griddata(lat,lon,GLA_,double(LAT),double(LON),'linear');
TUN_int    = griddata(lat,lon,TUN_,double(LAT),double(LON),'linear');
ROC_int    = griddata(lat,lon,ROC_,double(LAT),double(LON),'linear');
WET_int    = griddata(lat,lon,WET_,double(LAT),double(LON),'linear');
LAK_int    = griddata(lat,lon,LAK_,double(LAT),double(LON),'linear');
RIV_int    = griddata(lat,lon,RIV_,double(LAT),double(LON),'linear');
CLASS_int  = griddata(lat,lon,CLASS_,double(LAT),double(LON),'nearest','extrap', 'none'); % discrete value

land_mask  = NaN(size(GPP_TrendStore,1),size(GPP_TrendStore,2));
land_mask(~isnan(squeeze(GPP_TrendStore(:,:,1))))=1;
CLASS_int(isnan(land_mask))= NaN;
BOR_int(isnan(land_mask))= NaN;
GLA_int(isnan(land_mask))= NaN;
TUN_int(isnan(land_mask))= NaN;
ROC_int(isnan(land_mask))= NaN;
WET_int(isnan(land_mask))= NaN;
LAK_int(isnan(land_mask))= NaN;
RIV_int(isnan(land_mask))= NaN;


% Save in a NetCDF file
if save_flag
    % Create the NetCDF file
    filename = [path_in,'BAWLD_dataset.nc'];
    ncid = netcdf.create(filename, 'NETCDF4');
    
    % Define the dimensions
    dim_1 = netcdf.defDim(ncid,'dim_1', size(LAT, 1));
    dim_2 = netcdf.defDim(ncid,'dim_2', size(LAT, 2));
    
    bor_id = netcdf.defVar(ncid,'BOR_int', 'double', [dim_1, dim_2]);
    netcdf.putAtt(ncid, bor_id, 'description', 'Percent of boreal forest within grid cell (%)');
    gla_id = netcdf.defVar(ncid,'GLA_int', 'double', [dim_1, dim_2]);
    netcdf.putAtt(ncid, gla_id, 'description', 'Percent of galciers and permanent snow within grid cell');
    tun_id = netcdf.defVar(ncid,'TUN_int', 'double', [dim_1, dim_2]);
    netcdf.putAtt(ncid, tun_id, 'description', 'Percent of tundra within grid cell');
    wet_id = netcdf.defVar(ncid,'WET_int', 'double', [dim_1, dim_2]);
    netcdf.putAtt(ncid, wet_id, 'description', 'Percent of wetland within grid cell');
    lak_id = netcdf.defVar(ncid,'LAK_int', 'double', [dim_1, dim_2]);
    netcdf.putAtt(ncid, lak_id, 'description', 'Percent of lakes within grid cell');
    riv_id = netcdf.defVar(ncid,'RIV_int', 'double', [dim_1, dim_2]);
    netcdf.putAtt(ncid, riv_id, 'description', 'Percent of rivers within grid cell');
    roc_id = netcdf.defVar(ncid,'ROC_int', 'double', [dim_1, dim_2]);
    netcdf.putAtt(ncid, roc_id, 'description', 'Percent of rockland within grid cell');
    class_id = netcdf.defVar(ncid,'CLASS_int', 'double', [dim_1, dim_2]);
    netcdf.putAtt(ncid, class_id, 'description', 'k-means classification');
    lat_id = netcdf.defVar(ncid,'LAT', 'double', [dim_1, dim_2]);
    netcdf.putAtt(ncid, lat_id, 'description', 'latitude grid cell');
    lon_id =netcdf.defVar(ncid,'LON', 'double', [dim_1, dim_2]);
    netcdf.putAtt(ncid, lon_id, 'description', 'longitude grid cell');

    
    netcdf.putVar(ncid,bor_id,BOR_int);
    netcdf.putVar(ncid,gla_id,GLA_int);
    netcdf.putVar(ncid,tun_id,TUN_int);
    netcdf.putVar(ncid,wet_id,WET_int);
    netcdf.putVar(ncid,lak_id,LAK_int);
    netcdf.putVar(ncid,riv_id,RIV_int);
    netcdf.putVar(ncid,roc_id,ROC_int);
    netcdf.putVar(ncid,class_id,CLASS_int);
    netcdf.putVar(ncid,lat_id,LAT); 
    netcdf.putVar(ncid,lon_id,LON); 
    
    % Close the NetCDF file
    netcdf.close(ncid);
    
end

%% ------------------------------------------------
% -------------------------------------------------
% PLOTTING routines
% -------------------------------------------------
% -------------------------------------------------

% Example of boreal forest percent 
f = figure;
f.Position = [100 100 1700 1600];
title_list = {'galciers and permanent snow','rockland','tundra',...
              'boreal forest','wetland','lakes','rivers '};
          
elem_list  = {'GLA_int','ROC_int','TUN_int','BOR_int','WET_int','LAK_int','RIV_int'}; 
for i=1:7
    subplot(2,4,i)
    m_proj('Azimuthal Equal-area', 'lat',90,'long',30,'radius',30);
    eval(['m_pcolor(double(LON), double(LAT),',elem_list{i},');'])
    shading flat;
    c = colorbar;
    % Adjust colorbar position
    set(c, 'Location', 'southoutside','Fontsize',18)
    ylabel(c,'Percent (%)')
    m_coast('color', 'k');
    m_grid('box', 'fancy', 'fontsize', 16);
    t = title(title_list{i});  
    set(t,'Position',get(t,'Position')+[0 .03 0],'Fontsize',18);  % move up slightly
    set(gca,'Fontsize',18)
end

if save_flag
    saveas(gcf,[path_fig,'BAWLD_percents_classification.png'])
    saveas(gcf,[path_fig,'BAWLD_percents_classification.fig'])
end


%% -------------------------------------------------
% Create a figure with the classification of BAWLD
% -------------------------------------------------
f = figure;
f.Position = [100 100 1700 600];
% Plot the map of the variable you are interested in (e.g., 'Population')
m_proj('Azimuthal Equal-area', 'lat',90,'long',30,'radius',30); 
m_pcolor(double(LON), double(LAT), CLASS_int); 
colorbar;

n = 15; % Number of colors
% Define custom colormap with high contrast
colors = [
    0.1216, 0.4667, 0.7059;  % Blue Permafrost Peatlands
    1.0000, 0.4980, 0.0549;  % Orange Sparse Boreal Peatlands
    0.1725, 0.6275, 0.1725;  % Green Rivers
    0.8392, 0.1529, 0.1569;  % Red Glaciers
    0.5804, 0.4039, 0.7412;  % Purple Upland Tundra
    0.5490, 0.3373, 0.2941;  % Brown Common Boreal Peatlands
    0.8902, 0.4667, 0.7608;  % Pink Large Lakes
    0.4980, 0.4980, 0.4980;  % Gray Lake-rich wetlands
    0.7373, 0.7412, 0.1333;  % Yellow Dominant Boreal Peatlands
    0.0902, 0.7451, 0.8118;  % Cyan Wetland-rich Tundra
    0.9490, 0.3490, 0.2000;  % Dark Orange Alpine and Tundra Barrens
    0.9961, 0.5137, 0.7333;  % Light Pink Wetland and Lake-rich Tundra
    0.6078, 0.7490, 0.3569;  % Light Green Lake-rich Shield
    0.1804, 0.3490, 0.6706;  % Deep Blue Upland Boreal
    0.5686, 0.2353, 0.0980   % Dark Brown Wetland and Lake-rich Yedoma Tundra
];

% Normalize the colormap to range [0, 1]
colors = colors(1:n, :);
colors = colors ./ max(colors, [], 1);
colorbar;
colormap(colors)
c = colorbar;
caxis([1 n+1]); % Set the color range

% Adjust colorbar position
cb_pos    = get(c, 'Position');
cb_pos(1) = cb_pos(1) + 0; % Adjust the amount of shift as needed
set(c, 'Position', cb_pos,'Fontsize',18)

% Set the ticks to be at the centers of the color range
c.Ticks = 1 + (0.5:(n-0.5));

% Define the tick labels
% tickLabels = cellstr(num2str((1:n)', '%d'));
tickLabels = {'Permafrost Peatlands','Sparse Boreal Peatlands','Rivers',...
              'Glaciers','Upland Tundra','Common Boreal Peatlands','Large Lakes',...
              'Lake-rich wetlands','Dominant Boreal Peatlands','Wetland-rich Tundra',...
              'Alpine and Tundra Barrens','Wetland and Lake-rich Tundra','Lake-rich Shield','Upland Boreal','Wetland and Lake-rich Yedoma Tundra'};

% Assign the tick labels to the colorbar
c.TickLabels = tickLabels;

% Add a title to the colorbar
% c.Label.String = 'Discrete Colorbar';

% Add a title to the figure
t = title('BAWLD  classification');
set(t,'Position',get(t,'Position')+[0 .03 0],'Fontsize',18);  % move up slightly

% Optional: Adjust figure properties if needed
set(gca,'Fontsize',18)

m_coast('LineStyle','-','color','k');
m_grid('box','fancy','tickdir','in','fontsize', 18);

if save_flag
    saveas(gcf,[path_fig,'BAWLD_classification.png'])
    saveas(gcf,[path_fig,'BAWLD_classification.fig'])
end



% ----------------------------------------------------
% These lines compute the area distribution. It is not used for
% the classification map. 
% ----------------------------------------------------

% % Define the grid resolution (0.5 degrees in this case)
% resolution = 0.5;
% 
% % Create a uniform grid of latitudes and longitudes
% latitudes = 50:resolution:90;
% longitudes = -180:resolution:180;
% 
% % Define the projection as North Pole Lambert Azimuthal Equal Area
% projection = 'lambertazimuthal';
% latLonOrigin = [90 0]; % North Pole
% 
% % Define the radius of the Earth (WGS84 ellipsoid semi-major axis)
% radius = 6378137; % in meters
% 
% % Calculate the central latitude and longitude for each cell
% [gridLon, gridLat] = meshgrid(longitudes, latitudes);
% 
% % Preallocate an array to store the values for each grid cell
% gridValues = zeros(size(gridLon));
% 
% % Loop through each polygon in the shapefile
% for i = 1:length(S)
%   % Get the vertices of the polygon in meters
%     polyX = S(i).X;
%     polyY = S(i).Y;
%     
%     % Convert the polygon vertices to latitudes and longitudes
%     polyLon = zeros(size(polyX));
%     polyLat = zeros(size(polyY));
%     for j = 1:numel(polyX)
%         [polyLat(j), polyLon(j)] = projfwd_custom(projection, latLonOrigin, radius, polyY(j), polyX(j));
%     end
%     
%     % Get the area of the polygon
%     polyArea = S(i).Shp_Area;
%     
%     % Check if the polygon is empty (has no vertices)
%     if isempty(polyLon) || isempty(polyLat)
%         continue;
%     end
%     
%     % Check if the polygon is not a hole
%     if ~isnan(polyLon(1)) && ~isnan(polyLat(1))
%         % Find grid cells that contain the polygon
%         inPolygon = inpolygon(gridLon, gridLat, polyLon, polyLat);
%         
%         % Add the area of the polygon to the corresponding grid cells
%         gridValues(inPolygon) = gridValues(inPolygon) + polyArea;
%     end
% end
% 
% % Display or save the resulting grid
% % For example, to display the grid as an image:
% figure;
% imagesc(longitudes, latitudes, gridValues);
% colormap jet;
% colorbar;
% xlabel('Longitude');
% ylabel('Latitude');
% title('Grid with Polygon Areas');