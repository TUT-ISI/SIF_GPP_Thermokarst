% ###################################################################
% ####### COMPUTE ANOMALIES CSIF DATA  ####
% This script reads the cSIF database generated with OCO-2 SIF data and 
% MODIS surface reflectance every 4-day at 0.05 degrees. 
% This database also disentangles the sunlit and all sky SIF prediction
% 
% %Zhang, Y., Joiner, J., Alemohammad, S. H., Zhou, S., & Gentine, P. (2018).
% A global spatially contiguous solar-induced fluorescence (CSIF) dataset 
% using neural networks. Biogeosciences, 15(19), 5779-5800.
% https://bg.copernicus.org/articles/15/5779/2018/

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


path_in    = config.input_dir.cSIF; 
path_inHR  = config.input_dir.cSIFHR; 
path_ingeo = config.input_dir.GEO; 

save_flag   = 1;
years_      = 'cSIF_clearHR';
path_fig    = config.output_dir.general; 


% read the latitude and longitude range of interest for the Arctic Area
latitude  = ncread([path_ingeo,'latlon_Arctic.nc'],'latitude');
longitude = ncread([path_ingeo,'latlon_Arctic.nc'],'longitude');
[LAT,LON] = meshgrid(latitude,longitude); % coordinates of the SIF and GPP data 


switch years_
    case 'initial'
        year_ini = 2002; year_end = 2016;
        num_y    = year_end-year_ini+1;
    case 'final'
        year_ini = 2003; year_end = 2018;
        num_y    = year_end-year_ini+1;
    case 'cSIF_all'
        year_ini = 2003; year_end = 2016;
        num_y    = year_end-year_ini+1;
        
    case {'cSIF_clear','cSIF_clearHR'}
        year_ini = 2003; year_end = 2018;
        num_y    = year_end-year_ini+1;
end


% ncdisp([path_in,'OCO2.SIF.all.daily.',num2str(year_i),'.nc'])
% year_i = 2001;

% summer evaluation
% 1 June : 31 July 

for year_i = year_ini:year_end
    year_i
    
    switch years_
        
        case 'cSIF_all'
            lat = ncread([path_in,'OCO2.SIF.all.daily.',num2str(year_i),'.nc'],'lat');
            lon = ncread([path_in,'OCO2.SIF.all.daily.',num2str(year_i),'.nc'],'lon');
            [LAT_cSIF,LON_cSIF] = meshgrid(lat,lon); % coordinates of the cSIF data
            LAT_cSIF = double(LAT_cSIF);
            LON_cSIF = double(LON_cSIF);
            doy = ncread([path_in,'OCO2.SIF.all.daily.',num2str(year_i),'.nc'],'doy');
            all_sif = ncread([path_in,'OCO2.SIF.all.daily.',num2str(year_i),'.nc'],'all_daily_sif');
            
        case 'cSIF_clear'
            lat = ncread([path_in,'OCO2.SIF.clear.daily.',num2str(year_i),'.nc'],'lat');
            lon = ncread([path_in,'OCO2.SIF.clear.daily.',num2str(year_i),'.nc'],'lon');
            [LAT_cSIF,LON_cSIF] = meshgrid(lat,lon); % coordinates of the cSIF data
            LAT_cSIF = double(LAT_cSIF);
            LON_cSIF = double(LON_cSIF);
            doy = ncread([path_in,'OCO2.SIF.clear.daily.',num2str(year_i),'.nc'],'doy');
            all_sif = ncread([path_in,'OCO2.SIF.clear.daily.',num2str(year_i),'.nc'],'clear_daily_sif');
            
        case 'cSIF_clearHR' % Higher resolution 0.05 degrees   
            d_1 =  datetime(year_i,06,01); doy_1 = day(d_1,'dayofyear');
            d_2 =  datetime(year_i,07,31); doy_2 = day(d_2,'dayofyear');
            doy =  1:4:365;
            [~,pos] = find(doy <= doy_2 & doy >= doy_1);           
            
            for j = 1:length(pos)                
                if j==1
                    lat = ncread([path_inHR,num2str(year_i),filesep,'OCO2.SIF.clear.inst.',num2str(year_i),num2str(doy(pos(j))),'.v2.nc'],'lat');
                    lon = ncread([path_inHR,num2str(year_i),filesep,'OCO2.SIF.clear.inst.',num2str(year_i),num2str(doy(pos(j))),'.v2.nc'],'lon');
                    [LAT_cSIF,LON_cSIF] = meshgrid(lat,lon); % coordinates of the cSIF data HR
                    LAT_cSIF = double(LAT_cSIF);
                    LON_cSIF = double(LON_cSIF);
                end
                all_sif = ncread([path_inHR,num2str(year_i),filesep,'OCO2.SIF.clear.inst.',num2str(year_i),num2str(doy(pos(j))),'.v2.nc'],'clear_daily_SIF'); % There is also clear_inst_SIF
                csif_summer(:,:,j) = griddata(LAT_cSIF,LON_cSIF,double(all_sif),LAT,LON,'linear');
            end
            eval(['cSIF_SUM.year_',num2str(year_i),' = csif_summer;'])
            
            
    end
    
    switch years_
        case {'cSIF_all', 'cSIF_clear'}
            
            d_1 =  datetime(year_i,06,01); doy_1 = day(d_1,'dayofyear');
            d_2 =  datetime(year_i,07,31); doy_2 = day(d_2,'dayofyear');
            
            [pos,~] = find(doy <= doy_2 & doy >= doy_1);
            
            sif_aux     = double(all_sif(:,:,pos));
            csif_summer = NaN(size(LAT,1),size(LAT,2),size(sif_aux,3));
            for j = 1:size(sif_aux,3)
                csif_summer(:,:,j) = griddata(LAT_cSIF,LON_cSIF,squeeze(sif_aux(:,:,j)),LAT,LON,'linear');
            end
            eval(['cSIF_SUM.year_',num2str(year_i),' = csif_summer;'])
            clear csif_summer            
    end
    
end 

switch years_
    case 'cSIF_all'
        cSIF_SUM_ALL = cat(3,cSIF_SUM.year_2003,cSIF_SUM.year_2004,cSIF_SUM.year_2005,...
            cSIF_SUM.year_2006,cSIF_SUM.year_2007,cSIF_SUM.year_2008,cSIF_SUM.year_2009,...
            cSIF_SUM.year_2010,cSIF_SUM.year_2011,cSIF_SUM.year_2012,cSIF_SUM.year_2013,...
            cSIF_SUM.year_2014,cSIF_SUM.year_2015,cSIF_SUM.year_2016);
    case {'cSIF_clear','cSIF_clearHR'}
        cSIF_SUM_ALL = cat(3,cSIF_SUM.year_2003,cSIF_SUM.year_2004,cSIF_SUM.year_2005,...
            cSIF_SUM.year_2006,cSIF_SUM.year_2007,cSIF_SUM.year_2008,cSIF_SUM.year_2009,...
            cSIF_SUM.year_2010,cSIF_SUM.year_2011,cSIF_SUM.year_2012,cSIF_SUM.year_2013,...
            cSIF_SUM.year_2014,cSIF_SUM.year_2015,cSIF_SUM.year_2016,cSIF_SUM.year_2017,...
            cSIF_SUM.year_2018);        
end

mean_cSIF_sum   = mean(cSIF_SUM_ALL,3,'omitnan');
median_cSIF_sum = median(cSIF_SUM_ALL,3,'omitnan');

cSIF_mean   = mean_cSIF_sum;
cSIF_median = median_cSIF_sum;

cont = 0;
for year_i = year_ini:year_end
    cont = cont+1;
    eval(['ano_mean_sif   = mean(cSIF_SUM.year_',num2str(year_i),'-mean_cSIF_sum,3,''omitnan'');'])
    eval(['ano_median_sif = mean(cSIF_SUM.year_',num2str(year_i),'-median_cSIF_sum,3,''omitnan'');'])
    cSIF_anomaly_mean(cont,:,:)   = ano_mean_sif;
    cSIF_anomaly_median(cont,:,:) = ano_median_sif;
end


switch years_
    case 'cSIF_all'
        years_covered_cSIF = [2003:2016];        
        save([path_ingeo,'cSIFAnomaly_JunJul_05deg_2003_2016.mat'],'cSIF_anomaly_mean','cSIF_anomaly_median',...
            'cSIF_median','cSIF_mean','longitude','latitude','years_covered_cSIF')
        
    case 'cSIF_clear'
        years_covered_cSIF = [2003:2018];
        save([path_ingeo,'cSIFclearAnomaly_JunJul_05deg_2003_2018.mat'],'cSIF_anomaly_mean','cSIF_anomaly_median',...
            'cSIF_median','cSIF_mean','longitude','latitude','years_covered_cSIF')
        
    case 'cSIF_clearHR'
        years_covered_cSIF = [2003:2018];
        save([path_ingeo,'cSIFclearAnomaly_JunJul_005deg_2003_2018.mat'],'cSIF_anomaly_mean','cSIF_anomaly_median',...
            'cSIF_median','cSIF_mean','longitude','latitude','years_covered_cSIF')
end











