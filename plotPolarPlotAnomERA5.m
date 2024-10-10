% ##########################################################
% ####### PLOT POLAT PLOTS ANOMALIES ERA5   ###############
% This script generates polar plots for the 
% ECMWF 8-day interval data used in the 'SHAP 
% analyses' along different years 
% The ERA5 data has the following climate-realted variables:


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
path_in     = config.input_dir.ERA5_Land8day;

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
filename    = ['aggregateddata_',num2str(years_num(1)),'.nc'];
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
slhf        = NaN(num_years,longitude,latitude,8); % 8-day interval for 2 months (Jun-Jul) means max 8 observations
ssr         = slhf;
sshf        = slhf;
ssrd        = slhf;
tp          = slhf;
sp          = slhf;
lai_lv      = slhf;
lai_hv      = slhf;
fal         = slhf;

slhf_anom    = NaN(num_years,longitude,latitude);
ssr_anom     = slhf_anom;
sshf_anom    = slhf_anom;
tp_anom      = slhf_anom;
sp_anom      = slhf_anom;
lai_lv_anom  = slhf_anom;
lai_hv_anom  = slhf_anom;
fal_anom     = slhf_anom;

var_list_names       = {'fal','lai_hv','lai_lv','slhf','ssr','sp','sshf','ssrd','tp'};
var_list_names_title = {'fal','lai\_hv','lai\_lv','slhf','ssr','sp','sshf','ssrd','tp'};
units_list           = {'[-]','[-]','[-]','J.m^{-2}','J.m^{-2}','Pa','J.m^{-2}','J.m^{-2}','m'};
var_list             = length(var_list_names);

for i = 1:num_years
    filename                    = ['aggregateddata_',num2str(years_num(i)),'.nc'];
    
    time                        = ncread([path_in,filename],'time'); % 'days since 20XX-01-01 00:00:00'
    dateStr                     = ['1 June ',num2str(years_num(i))];    % Initial days of the year for June and July
    dateObj                     = datetime(dateStr, 'InputFormat', 'd MMMM yyyy');
    dayOfYear_start             = day(dateObj, 'dayofyear');    
    dateStr                     = ['31 July ',num2str(years_num(i))];   % Final days of the year for June and July
    dateObj                     = datetime(dateStr, 'InputFormat', 'd MMMM yyyy');
    dayOfYear_end               = day(dateObj, 'dayofyear');
    [pos,~]                     = find(time>=dayOfYear_start & time <=dayOfYear_end); % Values in June and July 
     
    % I have to check if there is something wrong in the netCDF aggregated
    fal(i,:,:,1:length(pos))    = ncread([path_in,filename],'fal',[1,1,pos(1)],[longitude, latitude,length(pos)]);
    lai_hv(i,:,:,1:length(pos)) = ncread([path_in,filename],'lai_hv',[1,1,pos(1)],[longitude, latitude,length(pos)]);
    lai_lv(i,:,:,1:length(pos)) = ncread([path_in,filename],'lai_lv',[1,1,pos(1)],[longitude, latitude,length(pos)]);
    slhf(i,:,:,1:length(pos))   = ncread([path_in,filename],'slhf',[1,1,pos(1)],[longitude, latitude,length(pos)]); % 'J m**-2' 'surface_upward_latent_heat_flux'
    ssr(i,:,:,1:length(pos))    = ncread([path_in,filename],'ssr',[1,1,pos(1)],[longitude, latitude,length(pos)]);
    sp(i,:,:,1:length(pos))     = ncread([path_in,filename],'sp',[1,1,pos(1)],[longitude, latitude,length(pos)]);
    sshf(i,:,:,1:length(pos))   = ncread([path_in,filename],'sshf',[1,1,pos(1)],[longitude, latitude,length(pos)]);
    ssrd(i,:,:,1:length(pos))   = ncread([path_in,filename],'ssrd',[1,1,pos(1)],[longitude, latitude,length(pos)]);
    tp(i,:,:,1:length(pos))     = ncread([path_in,filename],'tp',[1,1,pos(1)],[longitude, latitude,length(pos)]);

end 


%% Compute the mean for all the years in Jun-Jul
fal      = permute(fal,[1,4,2,3]);
lai_hv   = permute(lai_hv,[1,4,2,3]);
lai_lv   = permute(lai_lv,[1,4,2,3]);
slhf     = permute(slhf,[1,4,2,3]);
ssr      = permute(ssr,[1,4,2,3]);
sp       = permute(sp,[1,4,2,3]);
sshf     = permute(sshf,[1,4,2,3]);
ssrd     = permute(ssrd,[1,4,2,3]);
tp       = permute(tp,[1,4,2,3]);

for k =1:var_list 
    eval([var_list_names{k},'_mean = mean(reshape(',var_list_names{k},',[size(',var_list_names{k},',1).*size(',var_list_names{k},',2),size(',var_list_names{k},',3),size(',var_list_names{k},',4)]),1,''omitnan'');'])
end 


for i=1:num_years
    for k = 1:var_list 
    eval(['aux = squeeze(',var_list_names{k},'(i,:,:,:))-',var_list_names{k},'_mean;']) 
    eval([var_list_names{k},'_anom(i,:,:) = mean(aux,1,''omitnan'');'])
    end
end

% Set the limits for the figure
lim_1 = [-0.3,-0.07,-15e-3,-2e6,-2e6,-600,-2e6,-3.5e6,-2e-3]; %'fal','lai_hv','lai_lv','slhf','ssr','sp','sshf','ssrd','tp'};
lim_2 = [0.3,0.07,0.06,2e6,2e6,600,2e6,3.5e6,2e-3];


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
            caxis([lim_1(k),lim_2(k)])
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
                    saveas(gcf,[path_fig,'ERA5_anomalies',filesep,var_list_names{k},'_anomaly_Arctic_005deg_2002_2016_lim.png'])
                    saveas(gcf,[path_fig,'ERA5_anomalies',filesep,var_list_names{k},'_anomaly_Arctic_005deg_2002_2016_lim.fig'])
                case 'final'
                    saveas(gcf,[path_fig,'ERA5_anomalies',filesep,var_list_names{k},'_anomaly_Arctic_005deg_2003_2018_lim.png'])
                    saveas(gcf,[path_fig,'ERA5_anomalies',filesep,var_list_names{k},'_anomaly_Arctic_005deg_2003_2018_lim.fig'])
                    years_cov = [2003:2016];
                    
                case 'cSIF_all'
                    saveas(gcf,[path_fig,'ERA5_anomalies',filesep,var_list_names{k},'_anomaly_Arctic_005deg_2003_2016_lim.png'])
                    saveas(gcf,[path_fig,'ERA5_anomalies',filesep,var_list_names{k},'_anomaly_Arctic_005deg_2003_2016_lim.png'])
                    years_cov = [2003:2016];
                    
            end
        end
    end
    
end

%% Resample the data at the ArcticSIF resolution of 0.05
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
            caxis([lim_1(k),lim_2(k)])
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

% Save all the variables used in ERA5-Land and ARCLIM with the resolution of the
% ArcticSIF studies in a netCDF format
% NetCDF file
switch years_
    case 'initial'
        filename = [path_out,'ERA5_mean_anomalies_2002_2016.nc'];
    case 'final'
        filename = [path_out,'ERA5_mean_anomalies_2003_2018.nc'];
    case 'cSIF_all'
        filename = [path_out,'ERA5_mean_anomalies_2003_2016.nc'];
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
    eval(['netcdf.putAtt(ncid,varid_',var_list_names{k},'_mean_r,''description'',''ERA5 variable at 8-day interval mean for Jun-Jul for the years in the file name'');'])
    eval(['netcdf.putAtt(ncid,varid_',var_list_names{k},'_anom_r,''description'',''ERA5 variable at 8-day interval anomaly for Jun-Jul for each of the years in the file name'');'])
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




