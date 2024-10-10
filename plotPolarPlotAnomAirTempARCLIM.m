% ##########################################################
% ####### PLOT POLAT PLOTS ANOMALIES ARCLIM   ###############
% This script generates polar plots for the 
% ARCLIM 8-day interval data used in the 'SHAP 
% analyses' along different years 
% The ARCLIM data has the following climate-realted variables:
% - air Temperature
% - VPD

% Author: Neus Sabater
% e-mail: neus.sabater@fmi.fi 
% version 0: (12/06/2024) 
% ------------------------------------------------------
% ------------------------------------------------------
% setting for saving and ploting 

plotting    = 0;
save_flag   = 0;
years_      = 'final';  

addpath(genpath([pwd,filesep,'functions_environment']))
% Automatically detect environment
env = detect_environment();
% Load configuration
config = load_conf(env);


path_fig    = config.output_dir.general; 
path_out    = config.input_dir.computed; % The outputs are inputs for other scripts, so it goes to the computed directory

% ------------------------------------------------------
% Read data from ECMWF
% ------------------------------------------------------
% ------------------------------------------------------
path_in     = config.input_dir.ARCLIM;

switch years_
    case 'initial'
        years_num = [2002:2016];
        num_years = length(years_num);
        
    case 'final'
        years_num = [2003:2018];
        num_years = length(years_num);
        
    case 'cSIF_all'
        years_num = [2003:2016];
        num_years = length(years_num);
end

   
% Read dimensions (equal in all files)
filename    = ['2t_DMEA_',num2str(years_num(1)),'.nc'];
h           = ncinfo([path_in,filename]);
for j = 1:length(h.Dimensions)
    dimname = h.Dimensions(j).Name;
    dimlen  = h.Dimensions(j).Length;
    varname = matlab.lang.makeValidName(dimname);    % Ensure the variable name is valid by replacing invalid characters
    eval([varname ' = dimlen;']);    % Dynamically create the variable and assign its length
end

lat         = ncread([path_in,filename],'latitude');
lon         = ncread([path_in,filename],'longitude');

% Memory allocation
t2m       = NaN(num_years,longitude,latitude,8); % 8-day interval for 2 months (Jun-Jul) means max 8 observations
vpd       = t2m;

t2m_anom    = NaN(num_years,longitude,latitude);
vpd_anom    = t2m_anom;


for i = 1:num_years
    filename                    = ['2t_DMEA_',num2str(years_num(i)),'.nc'];
    filename_vpd                = ['vpd_DMEA_',num2str(years_num(i)),'.nc'];
       
    time                        = ncread([path_in,filename],'time'); % 'days since 20XX-01-01 00:00:00'
    dateStr                     = ['1 June ',num2str(years_num(i))];    % Initial days of the year for June and July
    dateObj                     = datetime(dateStr, 'InputFormat', 'd MMMM yyyy');
    dayOfYear_start             = day(dateObj, 'dayofyear');    
    dateStr                     = ['31 July ',num2str(years_num(i))];   % Final days of the year for June and July
    dateObj                     = datetime(dateStr, 'InputFormat', 'd MMMM yyyy');
    dayOfYear_end               = day(dateObj, 'dayofyear');
    [pos,~]                     = find(time>=dayOfYear_start & time <=dayOfYear_end); % Values in June and July 
     
    % I have to check if there is something wrong in the netCDF aggregated
    t2m(i,:,:,1:length(pos))    = ncread([path_in,filename],'t2m',[1,1,pos(1)],[longitude, latitude,length(pos)]);
    vpd(i,:,:,1:length(pos))    = ncread([path_in,filename_vpd],'vpd',[1,1,pos(1)],[longitude, latitude,length(pos)]);
    
end 

%% Compute the mean for all the years in Jun-Jul
t2m      = permute(t2m,[1,4,2,3]);
vpd      = permute(vpd,[1,4,2,3]);


t2m_mean = mean(reshape(t2m,[size(t2m,1).*size(t2m,2),size(t2m,3),size(t2m,4)]),1,'omitnan');
vpd_mean = mean(reshape(vpd,[size(vpd,1).*size(vpd,2),size(vpd,3),size(vpd,4)]),1,'omitnan');




for i=1:num_years
    aux = squeeze(t2m(i,:,:,:))-t2m_mean; 
    t2m_anom(i,:,:) = mean(aux,1,'omitnan');
    
    aux = squeeze(vpd(i,:,:,:))-vpd_mean; 
    vpd_anom(i,:,:) = mean(aux,1,'omitnan');    
end

% Set the limits for the figure
var_list_names       = {'t2m','vpd'};
var_list_names_title = var_list_names;
var_list             = length(var_list_names);
lim_1                = [-5,-0.07];
lim_2                = [5,0.07];


addpath(genpath(config.matlab_paths.mappingToolbox))
if plotting == 1
    for k=1:var_list
        k
        f = figure;
        f.Position = [100 100 2000 2000];
        for i=1:num_years
            subplot(4,4,i)
            m_proj('stereographic', 'lat',90,'long',30,'radius',30);
            eval(['m_pcolor(lon,lat,squeeze(',var_list_names{k},'_anom(i,:,:))'')'])
            shading flat; % Optional: removes grid lines in pcolor plots
            m_coast('LineStyle','-','color','k');
            m_grid('box','fancy','tickdir','in','fontsize', 18);
            h         = colorbar;
%             caxis([lim_1(k),lim_2(k)])
            % Adjust colorbar position
            cb_pos    = get(h, 'Position');
            cb_pos(1) = cb_pos(1) + 0.02; % Adjust the amount of shift as needed
            set(h, 'Position', cb_pos,'Fontsize',18)
            ylabel(h,[var_list_names_title{k},' anomaly ',units_list{k}],'Fontsize',18)
            cmap        = redblue(length(-0.1:0.02:0.1));
            [row,~]     = find(sum(cmap,2)==3);
            cmap(row,:) = [0.6,0.6,0.6];
            colormap(cmap);
            switch years_
                case 'initial'
                    title([num2str(i+2001)],'Fontsize',18) % Started in 2002
                case 'final'
                    title([num2str(i+2002)],'Fontsize',18) % Started in 2003
                case 'cSIF_all'
                    title([num2str(i+2002)],'Fontsize',18) % Started in 2003
            end
        end
        
        
        if save_flag
            switch years_
                case 'initial'
                    saveas(gcf,[path_fig,'ARCLIM_anomalies',filesep,var_list_names{k},'_anomaly_Arctic_005deg_2002_2016_lim.png'])
                    saveas(gcf,[path_fig,'ARCLIM_anomalies',filesep,var_list_names{k},'_anomaly_Arctic_005deg_2002_2016_lim.fig'])
                case 'final'
                    saveas(gcf,[path_fig,'ARCLIM_anomalies',filesep,var_list_names{k},'_anomaly_Arctic_005deg_2003_2018_lim.png'])
                    saveas(gcf,[path_fig,'ARCLIM_anomalies',filesep,var_list_names{k},'_anomaly_Arctic_005deg_2003_2018_lim.fig'])
                    years_cov = [2003:2016];
                    
                case 'cSIF_all'
                    saveas(gcf,[path_fig,'ARCLIM_anomalies',filesep,var_list_names{k},'_anomaly_Arctic_005deg_2003_2016_lim.png'])
                    saveas(gcf,[path_fig,'ARCLIM_anomalies',filesep,var_list_names{k},'_anomaly_Arctic_005deg_2003_2016_lim.png'])
                    years_cov = [2003:2016];
                    
            end
        end
    end
    
end

%Resample the data at the ArcticSIF resolution of 0.05
% read the latitude and longitude range of interest for the Arctic Area
path_ingeo = config.input_dir.computed;
latitude   = ncread([path_ingeo,'latlon_Arctic.nc'],'latitude');
longitude  = ncread([path_ingeo,'latlon_Arctic.nc'],'longitude');
[LAT,LON]  = meshgrid(latitude,longitude); % coordinates of the SIF and GPP data 

[LAT_era5,LON_era5] = meshgrid(lat,lon); % coordinates of the ERA5 data



for k=1:var_list
    eval([var_list_names{k},'_mean_r = griddata(LAT_era5,LON_era5,squeeze(',var_list_names{k},'_mean(1,:,:)),LAT,LON,''linear'');']) 
    for i = 1:num_years
        eval([var_list_names{k},'_anom_r(i,:,:) = griddata(LAT_era5,LON_era5,squeeze(',var_list_names{k},'_anom(i,:,:)),LAT,LON,''linear'');'])
    end
end


if plotting == 1
    for k=1:var_list
        k
        f = figure;
        f.Position = [100 100 2000 2000];
        for i=1:num_years
            subplot(4,4,i)
            m_proj('stereographic', 'lat',90,'long',30,'radius',30);
            eval(['m_pcolor(longitude,latitude,squeeze(',var_list_names{k},'_anom_r(i,:,:))'')'])
            shading flat; % Optional: removes grid lines in pcolor plots
            m_coast('LineStyle','-','color','k');
            m_grid('box','fancy','tickdir','in','fontsize', 18);
            h         = colorbar;
%             caxis([lim_1(k),lim_2(k)])
            % Adjust colorbar position
            cb_pos    = get(h, 'Position');
            cb_pos(1) = cb_pos(1) + 0.02; % Adjust the amount of shift as needed
            set(h, 'Position', cb_pos,'Fontsize',18)
            ylabel(h,[var_list_names_title{k},' anomaly ',units_list{k}],'Fontsize',18)
            cmap        = redblue(length(-0.1:0.02:0.1));
            [row,~]     = find(sum(cmap,2)==3);
            cmap(row,:) = [0.6,0.6,0.6];
            colormap(cmap);
            switch years_
                case 'initial'
                    title([num2str(i+2001)],'Fontsize',18) % Started in 2002
                case 'final'
                    title([num2str(i+2002)],'Fontsize',18) % Started in 2003
                case 'cSIF_all'
                    title([num2str(i+2002)],'Fontsize',18) % Started in 2003
            end
        end
        
    end
end 


%% Save all the variables used in ERA5-Land and ARCLIM with the resolution of the
% ArcticSIF studies in a netCDF format
% NetCDF file
switch years_
    case 'initial'
        filename = [path_out,'ARCLIM_mean_anomalies_2002_2016.nc'];
    case 'final'
        filename = [path_out,'ARCLIM_mean_anomalies_2003_2018.nc'];
    case 'cSIF_all'
        filename = [path_out,'ARCLIM_mean_anomalies_2003_2016.nc'];
end

% Create a new NetCDF file format 4 to allocate more memory
ncid = netcdf.create(filename, 'NETCDF4');

% Define dimensions
dimid_lon = netcdf.defDim(ncid, 'lon', length(longitude));
dimid_lat = netcdf.defDim(ncid, 'lat', length(latitude));
dimid_t   = netcdf.defDim(ncid, 'time_steps', num_years);

% Define variables

for k=1:var_list 
    eval(['varid_',var_list_names{k},'_mean_r = netcdf.defVar(ncid, ''',var_list_names{k},'__mean_r'', ''NC_FLOAT'', [dimid_lon,dimid_lat]);'])
    eval(['varid_',var_list_names{k},'_anom_r = netcdf.defVar(ncid, ''',var_list_names{k},'_anom_r'', ''NC_FLOAT'', [dimid_t,dimid_lon,dimid_lat]);'])
end 
varid_latitude          = netcdf.defVar(ncid, 'latitude', 'NC_FLOAT', dimid_lat);
varid_longitude         = netcdf.defVar(ncid, 'longitude', 'NC_FLOAT', dimid_lon);



% Add attributes for descriptions
for k=1:var_list 
    eval(['netcdf.putAtt(ncid,varid_',var_list_names{k},'_mean_r,''description'',''ARCLIM variable at 8-day interval mean for Jun-Jul for the years in the file name'');'])
    eval(['netcdf.putAtt(ncid,varid_',var_list_names{k},'_anom_r,''description'',''ARCLIM variable at 8-day interval anomaly for Jun-Jul for each of the years in the file name'');'])
end 


% End definitions and enter data mode
netcdf.endDef(ncid);


% Store the data into the NetCDF variables
for k=1:var_list 
    eval(['netcdf.putVar(ncid,varid_',var_list_names{k},'_mean_r,',var_list_names{k},'_mean_r);'])
    eval(['netcdf.putVar(ncid,varid_',var_list_names{k},'_anom_r,',var_list_names{k},'_anom_r);'])
end 

% Close the NetCDF file
netcdf.close(ncid);
disp(['Saved NetCDF file'])


%% Saving these variables in the previously created ERA-5 netCDF

switch years_
    case 'initial'
        filename = [path_out,'ERA5_mean_anomalies_2002_2016.nc'];
    case 'final'
        filename = [path_out,'ERA5_mean_anomalies_2003_2018.nc'];
    case 'cSIF_all'
        filename = [path_out,'ERA5_mean_anomalies_2003_2016.nc'];
end

% Open the existing NetCDF file in write mode
ncid = netcdf.open(filename, 'NC_WRITE');

% Enter define mode
netcdf.reDef(ncid);

% Define dimensions
% dimid_lon = netcdf.defDim(ncid, 'lon', length(longitude));
% dimid_lat = netcdf.defDim(ncid, 'lat', length(latitude));
% dimid_t   = netcdf.defDim(ncid, 'time_steps', num_years);

% Define variables

for k=1:var_list 
    eval(['varid_',var_list_names{k},'_mean_r = netcdf.defVar(ncid, ''',var_list_names{k},'__mean_r'', ''NC_FLOAT'', [dimid_lon,dimid_lat]);'])
    eval(['varid_',var_list_names{k},'_anom_r = netcdf.defVar(ncid, ''',var_list_names{k},'_anom_r'', ''NC_FLOAT'', [dimid_t,dimid_lon,dimid_lat]);'])
end 

% Add attributes for descriptions
for k=1:var_list 
    eval(['netcdf.putAtt(ncid,varid_',var_list_names{k},'_mean_r,''description'',''ARCLIM variable at 8-day interval mean for Jun-Jul for the years in the file name'');'])
    eval(['netcdf.putAtt(ncid,varid_',var_list_names{k},'_anom_r,''description'',''ARCLIM variable at 8-day interval anomaly for Jun-Jul for each of the years in the file name'');'])
end 


% End definitions and enter data mode
netcdf.endDef(ncid);


% Store the data into the NetCDF variables
for k=1:var_list 
    eval(['netcdf.putVar(ncid,varid_',var_list_names{k},'_mean_r,',var_list_names{k},'_mean_r);'])
    eval(['netcdf.putVar(ncid,varid_',var_list_names{k},'_anom_r,',var_list_names{k},'_anom_r);'])
end 

% Close the NetCDF file
netcdf.close(ncid);
disp(['Saved NetCDF file'])

% Change the name of the ERA file to also indicate that has the ARCLIM
% variables!!! e.g., ERA5_ARCLIM_mean_anomalies_2003_2018.nc




