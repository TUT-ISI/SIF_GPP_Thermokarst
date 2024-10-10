% ###################################################################
% ####### PLOT POLAT PLOTS ANOMALIES AIR TEMPERATURE ####
% This script generates polar plots for air temperature anomalies data from
% ECMWF 8-day interval along different years 
% %
% Author: Neus Sabater
% e-mail: neus.sabater@fmi.fi 
% version 0: (19/02/2024) 
% ------------------------------------------------------
% ------------------------------------------------------
% setting for saving and ploting 

plotting    = 0;
save_flag   = 0;
years_      = 'cSIF_all';  

addpath(genpath([pwd,filesep,'functions_environment']))
% Automatically detect environment
env = detect_environment();
% Load configuration
config = load_conf(env);


path_fig    = config.output_dir.general; 
path_out    = config.input_dir.computed; % The outputs are inputs for other scripts, so it goes to the computed directory

% ------------------------------------------------------
% Read the air Temperature data from ECMWF
% ------------------------------------------------------
% ------------------------------------------------------
path_in   = config.input_dir.ERA5T2m; 

names = {'T2m_2002_2006_8days','T2m_2007_2011_8days','T2m_2012_2016_8days',...
         'T2m_2017_2021_8days'};
     
% ncdisp([path_in,names{i},'.nc'])
T2m_aux = NaN(4,12,1440,721);
month_name     = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
year_ini_files = [2002,2007,2012,2017];
year_end_files = [2006,2011,2016,2021];
cont = 0;


for i = 1:4
    % 'days since 2002-01-01 00:00:00'
    %  years each file
    time = ncread([path_in,names{i},'.nc'],'time');
    lat  = ncread([path_in,names{i},'.nc'],'lat');
    lon  = ncread([path_in,names{i},'.nc'],'lon');
    t2m  = ncread([path_in,names{i},'.nc'],'t2m');
    
    
    t1 = datetime(year_ini_files(i),1,1,0,0,0);
    % Add all the rest of the time vector
    for time_step =1:length(time)
        t2(time_step) = t1 + days(time(time_step));
    end
    
    for month_num =6:7 % Months of June and July
        eval(['days_month_',num2str(month_num),'= t2(month(t2)==month_num);']);
        eval(['aux = days_month_',num2str(month_num),';'])
        for time_step =1:length(aux)
            [~,pos] = find(t2 == aux(time_step));
            POS(time_step) = pos;
        end
        eval(['POS_month',num2str(month_num),'= POS;'])
        clear aux POS pos
    end
    
    num_years = year_end_files(i)- year_ini_files(i)+1;
    
    % group the days corresponding to June and July each year
    [~,nn_month6]    = find(diff(POS_month6)>1);
    [~,nn_month7]    = find(diff(POS_month7)>1);
    
    end_nn_month6    = [nn_month6,length(POS_month6)];
    ini_nn_month6    = [1,nn_month6+1];
    
    end_nn_month7    = [nn_month7,length(POS_month7)];
    ini_nn_month7    = [1,nn_month7+1];
    
    for nn_ = 1:num_years
        cont = cont +1;
        t2m_june(cont,:,:,:) = squeeze(t2m(:,:,POS_month6(ini_nn_month6(nn_):end_nn_month6(nn_))));
        t2m_july(cont,:,:,:) = squeeze(t2m(:,:,POS_month7(ini_nn_month7(nn_):end_nn_month7(nn_))));
    end
end

t2m_june = permute(t2m_june,[1,4,2,3]);
t2m_july = permute(t2m_july,[1,4,2,3]);


%% Crop depending on the years to be accounted for in the time series

switch years_
    case 'initial' % 2002-2016
        t2m_june = t2m_june(1:15,:,:,:);
        t2m_july = t2m_july(1:15,:,:,:);
        
    case 'final' % 2003-2018
        t2m_june = t2m_june(2:17,:,:,:);
        t2m_july = t2m_july(2:17,:,:,:);
        
    case 'cSIF_all' % 2003-2016
        t2m_june = t2m_june(2:15,:,:,:);
        t2m_july = t2m_july(2:15,:,:,:);
end




%% Compute the temperature Anomaly in Jun-Jul
t2m_june_reshape = reshape(t2m_june,[size(t2m_june,1)*size(t2m_june,2),size(t2m_june,3),size(t2m_june,4)]);
t2m_july_reshape = reshape(t2m_july,[size(t2m_june,1)*size(t2m_june,2),size(t2m_june,3),size(t2m_june,4)]);
t2m_summer       = [t2m_june_reshape;t2m_july_reshape];
t2m_mean_summer  = squeeze(mean(t2m_summer,1,'omitnan'));

for i=1:size(t2m_june,1)
    t2m_summer_year_aux            = squeeze([t2m_june(i,:,:,:),t2m_july(i,:,:,:)]);
    t2m_summer_anomaly_year(i,:,:) = squeeze(mean(t2m_summer_year_aux,1))-t2m_mean_summer;
    
end 
    
    
%% Temperature anomalies

if plotting == 1
    f = figure;
    f.Position = [100 100 2000 2000];
    for i=1:size(t2m_june,1)
        subplot(4,4,i)
        m_proj('stereographic', 'lat',90,'long',30,'radius',30);
        m_pcolor(lon,lat,squeeze(t2m_summer_anomaly_year(i,:,:))'); % Ensure lon and lat match the grid of land_cover data
        shading flat; % Optional: removes grid lines in pcolor plots
        m_coast('LineStyle','-','color','k');
        m_grid('box','fancy','tickdir','in','fontsize', 18);
        h         = colorbar;
        caxis([-5,5])
        % Adjust colorbar position
        cb_pos    = get(h, 'Position');
        cb_pos(1) = cb_pos(1) + 0.02; % Adjust the amount of shift as needed
        set(h, 'Position', cb_pos,'Fontsize',18)
        ylabel(h,'air T Anomaly','Fontsize',18)
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
                saveas(gcf,[path_fig,'Tanomaly_Arctic_005deg_2002_2016.png'])
                saveas(gcf,[path_fig,'Tanomaly_Arctic_005deg_2002_2016.fig'])
                years_cov = [2002:2016];
                save('airTAnomaly_JunJul_025deg_2002_2016.mat','years_cov','t2m_summer_anomaly_year','lon','lat')
                
            case 'final'
                saveas(gcf,[path_fig,'Tanomaly_Arctic_005deg_2003_2018.png'])
                saveas(gcf,[path_fig,'Tanomaly_Arctic_005deg_2003_2018.fig'])
                years_cov = [2003:2016];
                save('airTAnomaly_JunJul_025deg_2003_2018.mat','years_cov','t2m_summer_anomaly_year','lon','lat')
                
            case 'cSIF_all'
                saveas(gcf,[path_fig,'Tanomaly_Arctic_005deg_2003_2016.png'])
                saveas(gcf,[path_fig,'Tanomaly_Arctic_005deg_2003_2016.fig'])
                years_cov = [2003:2016];
                save([path_out,'airTAnomaly_JunJul_025deg_2003_2016.mat'],'years_cov','t2m_summer_anomaly_year','lon','lat')
                
        end
    end
end


%% Temperature anomalies at 0.05 degrees resolution
% read the latitude and longitude range of interest for the Arctic Area
path_ingeo = config.input_dir.computed;
latitude   = ncread([path_ingeo,'latlon_Arctic.nc'],'latitude');
longitude  = ncread([path_ingeo,'latlon_Arctic.nc'],'longitude');
[LAT,LON]  = meshgrid(latitude,longitude); % coordinates of the SIF and GPP data 

[LAT_era5,LON_era5] = meshgrid(lat,lon); % coordinates of the ERA5 data

for i=1:size(t2m_summer_anomaly_year,1)
    t2m_anomaly_mean(i,:,:) = griddata(LAT_era5,LON_era5,squeeze(t2m_summer_anomaly_year(i,:,:)),LAT,LON,'linear');
end



if save_flag
    switch years_
        case 'initial'
            saveas(gcf,[path_fig,'Tanomaly_Arctic_005deg_2002_2016.png'])
            saveas(gcf,[path_fig,'Tanomaly_Arctic_005deg_2002_2016.fig'])
            years_cov = [2002:2016];
            save([path_out,'airTAnomaly_JunJul_005deg_2002_2016.mat'],'years_cov','t2m_anomaly_mean','longitude','latitude')
            
        case 'final'
            saveas(gcf,[path_fig,'Tanomaly_Arctic_005deg_2003_2018.png'])
            saveas(gcf,[path_fig,'Tanomaly_Arctic_005deg_2003_2018.fig'])
            years_cov = [2003:2016];
            save([path_out,'airTAnomaly_JunJul_005deg_2003_2018.mat'],'years_cov','t2m_anomaly_mean','longitude','latitude')
            
        case 'cSIF_all'
            saveas(gcf,[path_fig,'Tanomaly_Arctic_005deg_2003_2016.png'])
            saveas(gcf,[path_fig,'Tanomaly_Arctic_005deg_2003_2016.fig'])
            years_cov = [2003:2016];
            save([path_out,'airTAnomaly_JunJul_005deg_2003_2016.mat'],'years_cov','t2m_anomaly_mean','longitude','latitude')
            
    end
end


