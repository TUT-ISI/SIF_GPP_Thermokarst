% ###################################################################
% ################## OCO2 SIF classification  #######################
% This script reads the OCO_2 data from the NetCDF downloaded for the 
% years 2015, 2016, and 2017 to validate the the increasing SIF during 
% 2016 observed with the GOSIF is corroborated and not an artifact 
% due to the training with meteorological datasets. 

% OCO2 data download from:
% https://disc.gsfc.nasa.gov/datasets/OCO2_L2_Lite_SIF_10r/summary?keywords=oco2%20sif%20lit
% OCO-2 Level 2 bias-corrected solar-induced fluorescence and other select fields from the IMAP-DOAS algorithm aggregated as daily files, Retrospective processing V11r (OCO2_L2_Lite_SIF) at GES


% Author: Neus Sabater
% e-mail: neus.sabater@fmi.fi 
% version 0: (19/03/2024) 
%% ---------------------------
%  set parameters 
% ---------------------------

plotting    = 1;
save_flag   = 0;

addpath(genpath([pwd,filesep,'functions_environment']))
% Automatically detect environment
env = detect_environment();
% Load configuration
config = load_conf(env);


path_fig    = config.output_dir.general; 
path_in     = config.input_dir.OCO2; 
years_      = 'cSIF_clear';
listing     = dir(path_in);

% Create a datetime object for January 1st, 1990
start_date        = datetime(1990, 1, 1, 0, 0, 0);



% % Read the coordinates used in ArcticSIF
path_geo = config.input_dir.GEO;
lat = ncread([path_geo,'latlon_Arctic.nc'],'latitude');  % 0.05
lon = ncread([path_geo,'latlon_Arctic.nc'],'longitude'); % 0.05
[LAT,LON] = meshgrid(lat,lon);



cont_2015 = 1;
cont_2016 = 1;
cont_2017 = 1;


for k = 4:length(listing)
    k
    file = [path_in,listing(k).name];
    
    SIF_740nm  = ncread(file,'SIF_740nm');
    SIF_Uncertainty_740nm = ncread(file,'SIF_Uncertainty_740nm');
    Daily_SIF_740nm = ncread(file,'Daily_SIF_740nm');
    Quality_Flag    = ncread(file,'Quality_Flag'); % 0 = best (passes quality control + cloud fraction = 0.0); 1 = good (passes quality control); 2 = bad (failed quality control); -1 = not investigated'

    Daily_SIF_757nm = ncread(file,'Daily_SIF_757nm');
    Daily_SIF_771nm = ncread(file,'Daily_SIF_771nm');
    Delta_Time      = ncread(file,'Delta_Time'); % seconds since 1 Jan 1990 
    % Add the seconds to the start date
    result_date     = start_date + seconds(Delta_Time);
    % Extract year, month, and day
    result_year     = year(result_date);
    result_month    = month(result_date);
    result_day      = day(result_date);
    
    Longitude_Corners = ncread(file,'Longitude_Corners');
    Latitude_Corners  = ncread(file,'Latitude_Corners');
    Longitude         = ncread(file,'Longitude');
    Latitude          = ncread(file,'Latitude');
    [val,~]           = find(Latitude>60); % ArcticSIF region
    
    
    % Re-gridding 
    Quality_Flag     = Quality_Flag(val);
    
    data_values      = SIF_740nm(val);
    data_values      = data_values(Quality_Flag == 0 | 1);
    
    data_val_SIF_740nm_unc      = SIF_Uncertainty_740nm(val);
    data_val_SIF_740nm_unc      = data_val_SIF_740nm_unc(Quality_Flag == 0 | 1);

    data_val_Daily_SIF_740      = Daily_SIF_740nm(val);
    data_val_Daily_SIF_740      = data_val_Daily_SIF_740(Quality_Flag == 0 | 1);
    
    data_val_Daily_SIF_757      = Daily_SIF_757nm(val);
    data_val_Daily_SIF_757      = data_val_Daily_SIF_757(Quality_Flag == 0 | 1);
    
    data_val_Daily_SIF_771      = Daily_SIF_771nm(val);
    data_val_Daily_SIF_771      = data_val_Daily_SIF_771(Quality_Flag == 0 | 1);
    
    
    Longitude_values = Longitude(val); Longitude_values = Longitude_values(Quality_Flag == 0 | 1);
    Latitude_values  = Latitude(val); Latitude_values  = Latitude_values(Quality_Flag == 0 | 1);
    
   
    
    % Preallocate matrix to hold indices
    ind_val = zeros(size(data_values, 1), 2);
    
    for i = 1:size(data_values, 1)
        % Find the index of the closest latitude in the grid
        [~, lat_idx] = min(abs(lat - Latitude_values(i)));
        
        % Find the index of the closest longitude in the grid
        [~, lon_idx] = min(abs(lon - Longitude_values(i)));
        
        % Store the indices
        ind_val(i, :) = [lon_idx lat_idx,];
    end
    
    % Create a matrix to hold the colocated values
    grid_data            = NaN(size(LAT));
    grid_SIF_740nm       = NaN(size(LAT));
    grid_SIF_740nm_unc   = NaN(size(LAT));
    grid_Daily_SIF_740   = NaN(size(LAT));
    grid_Daily_SIF_757   = NaN(size(LAT));
    grid_Daily_SIF_771   = NaN(size(LAT));
    
    % Assign the data values to the grid
    for i = 1:size(data_values, 1)
        grid_SIF_740nm(ind_val(i, 1), ind_val(i, 2))       = data_values(i);
        grid_SIF_740nm_unc(ind_val(i, 1), ind_val(i, 2))   = data_val_SIF_740nm_unc(i);
        grid_Daily_SIF_740(ind_val(i, 1), ind_val(i, 2))   = data_val_Daily_SIF_740(i);
        grid_Daily_SIF_757(ind_val(i, 1), ind_val(i, 2))   = data_val_Daily_SIF_757(i);
        grid_Daily_SIF_771(ind_val(i, 1), ind_val(i, 2))   = data_val_Daily_SIF_771(i);
    end
    
    
    if all(result_year == 2015 & (result_month == 6 | 7))
        D2015_Daily_SIF_740nm(cont_2015,:,:) = grid_Daily_SIF_740;
        D2015_Daily_SIF_757nm(cont_2015,:,:) = grid_Daily_SIF_757;
        D2015_Daily_SIF_771nm(cont_2015,:,:) = grid_Daily_SIF_771;
        D2015_SIF_Uncertainty_740nm(cont_2015,:,:) = grid_SIF_740nm_unc;
        D2015_SIF_740nm(cont_2015,:,:)             = grid_SIF_740nm;
        
        cont_2015 = cont_2015 + 1;
        
    elseif all(result_year == 2016 & (result_month == 6 | 7))
        
        D2016_Daily_SIF_740nm(cont_2016,:,:) = grid_Daily_SIF_740;
        D2016_Daily_SIF_757nm(cont_2016,:,:) = grid_Daily_SIF_757;
        D2016_Daily_SIF_771nm(cont_2016,:,:) = grid_Daily_SIF_771;
        D2016_SIF_Uncertainty_740nm(cont_2016,:,:) = grid_SIF_740nm_unc;
        D2016_SIF_740nm(cont_2016,:,:) = grid_SIF_740nm;
        
        cont_2016 = cont_2016 + 1;
        
    elseif all(result_year == 2017 & (result_month == 6 | 7))
        
        D2017_Daily_SIF_740nm(cont_2017,:,:) = grid_Daily_SIF_740;
        D2017_Daily_SIF_757nm(cont_2017,:,:) = grid_Daily_SIF_757;
        D2017_Daily_SIF_771nm(cont_2017,:,:) = grid_Daily_SIF_771;
        D2017_SIF_Uncertainty_740nm(cont_2017,:,:) = grid_SIF_740nm_unc;
        D2017_SIF_740nm(cont_2017,:,:) = grid_SIF_740nm;
        
        cont_2017  = cont_2017 + 1;
    end

    
    clear grid_SIF_740nm  grid_SIF_740nm_unc  grid_Daily_SIF_740 grid_Daily_SIF_757nm grid_Daily_SIF_771nm
    
end


% Display the result
% disp(result_date);

% Now do the average if there is more than one value per grid 
D2015_SIF_740nm             = mean(D2015_SIF_740nm,1,'omitnan');
D2015_Daily_SIF_740nm       = mean(D2015_Daily_SIF_740nm,1,'omitnan');
D2015_Daily_SIF_757nm       = mean(D2015_Daily_SIF_757nm,1,'omitnan');
D2015_Daily_SIF_771nm       = mean(D2015_Daily_SIF_771nm,1,'omitnan');
D2015_SIF_Uncertainty_740nm = mean(D2015_SIF_Uncertainty_740nm,1,'omitnan');

D2016_SIF_740nm             = mean(D2016_SIF_740nm,1,'omitnan');
D2016_Daily_SIF_740nm       = mean(D2016_Daily_SIF_740nm,1,'omitnan');
D2016_Daily_SIF_757nm       = mean(D2016_Daily_SIF_757nm,1,'omitnan');
D2016_Daily_SIF_771nm       = mean(D2016_Daily_SIF_771nm,1,'omitnan');
D2016_SIF_Uncertainty_740nm = mean(D2016_SIF_Uncertainty_740nm,1,'omitnan');

D2017_SIF_740nm             = mean(D2017_SIF_740nm,1,'omitnan');
D2017_Daily_SIF_740nm       = mean(D2017_Daily_SIF_740nm,1,'omitnan');
D2017_Daily_SIF_757nm       = mean(D2017_Daily_SIF_757nm,1,'omitnan');
D2017_Daily_SIF_771nm       = mean(D2017_Daily_SIF_771nm,1,'omitnan');
D2017_SIF_Uncertainty_740nm = mean(D2017_SIF_Uncertainty_740nm,1,'omitnan');


%%
addpath(genpath('/Users/neus/Dropbox/DOCUMENTS_MyMac_Backup/Documents/MATLAB/m_map'))
f = figure;
f.Position = [100 100 800 600];
m_proj('stereographic', 'lat',90,'long',30,'radius',30);
m_pcolor(LON,LAT,squeeze(D2015_SIF_740nm(1,:,:))); % Ensure lon and lat match the grid of land_cover data
shading flat; % Optional: removes grid lines in pcolor plots
m_coast('LineStyle','-','color','k');
m_grid('box','fancy','tickdir','in','fontsize', 18);
h         = colorbar;
% Adjust colorbar position
caxis([-1,3])
cb_pos    = get(h, 'Position');
cb_pos(1) = cb_pos(1) + 0.05; % Adjust the amount of shift as needed
set(h, 'Position', cb_pos,'Fontsize',18)
ylabel(h,'SIF @ 740 nm','Fontsize',18)
title('SIF OCO-2 Lite v11 at 740nm  2015 (June - July)')
set(gca,'Fontsize',18)
if flag_save
    saveas(gcf,[path_fig,'SIF_OCO2_Litev11_Arctic_005deg2015.png'])
    saveas(gcf,[path_fig,'SIF_OCO2_Litev11_Arctic_005deg2015.fig'])
end

f = figure;
f.Position = [100 100 800 600];
m_proj('stereographic', 'lat',90,'long',30,'radius',30);
m_pcolor(LON,LAT,squeeze(D2016_SIF_740nm(1,:,:))); % Ensure lon and lat match the grid of land_cover data
shading flat; % Optional: removes grid lines in pcolor plots
m_coast('LineStyle','-','color','k');
m_grid('box','fancy','tickdir','in','fontsize', 18);
h         = colorbar;
% Adjust colorbar position
caxis([-1,3])
cb_pos    = get(h, 'Position');
cb_pos(1) = cb_pos(1) + 0.05; % Adjust the amount of shift as needed
set(h, 'Position', cb_pos,'Fontsize',18)
ylabel(h,'SIF @ 740 nm','Fontsize',18)
title('SIF OCO-2 Lite v11 at 740nm  2016 (June - July)')
set(gca,'Fontsize',18)
if flag_save
    saveas(gcf,[path_fig,'SIF_OCO2_Litev11_Arctic_005deg2016.png'])
    saveas(gcf,[path_fig,'SIF_OCO2_Litev11_Arctic_005deg2016.fig'])
end



f = figure;
f.Position = [100 100 800 600];
m_proj('stereographic', 'lat',90,'long',30,'radius',30);
m_pcolor(LON,LAT,squeeze(D2017_SIF_740nm(1,:,:))); % Ensure lon and lat match the grid of land_cover data
shading flat; % Optional: removes grid lines in pcolor plots
m_coast('LineStyle','-','color','k');
m_grid('box','fancy','tickdir','in','fontsize', 18);
h         = colorbar;
% Adjust colorbar position
caxis([-1,3])
cb_pos    = get(h, 'Position');
cb_pos(1) = cb_pos(1) + 0.05; % Adjust the amount of shift as needed
set(h, 'Position', cb_pos,'Fontsize',18)
ylabel(h,'SIF @ 740 nm','Fontsize',18)
title('SIF OCO-2 Lite v11 at 740nm  2017 (June - July)')
set(gca,'Fontsize',18)
if flag_save
    saveas(gcf,[path_fig,'SIF_OCO2_Litev11_Arctic_005deg2017.png'])
    saveas(gcf,[path_fig,'SIF_OCO2_Litev11_Arctic_005deg2017.fig'])
end


% ------------------------------------------------------------------
% Histograms for the Siberia anomaly (60 degree E and 120 degree E)
% ------------------------------------------------------------------
figure % 740 nm
aux_1 = squeeze(D2015_SIF_740nm(1,4800:6000,:));
aux_2 = squeeze(D2016_SIF_740nm(1,4800:6000,:));
aux_3 = squeeze(D2017_SIF_740nm(1,4800:6000,:));
histogram(aux_1(:),'Normalization','pdf','FaceAlpha',0.3,'FaceColor','r')
hold on 
histogram(aux_2(:),'Normalization','pdf','FaceAlpha',0.3,'FaceColor','k')
histogram(aux_3(:),'Normalization','pdf','FaceAlpha',0.3,'FaceColor','g')
legend(['SIF 2015 - mean value:',num2str(round(mean(aux_1(:),'omitnan'),3))],...
       ['SIF 2016 - mean value:',num2str(round(mean(aux_2(:),'omitnan'),3))],...
       ['SIF 2017 - mean value:',num2str(round(mean(aux_3(:),'omitnan'),3))])
set(gca,'Fontsize',18)
grid on 
box on 
xlabel('SIF at 740 nm (mW/(m^2 sr nm))')
ylabel('normalized distribution (pdf)')
title('SIF Area: 60^{\circ}N - 90^{\circ}N, 60^{\circ}E - 120^{\circ}E, ')
if flag_save
    saveas(gcf,[path_fig,'SIF_OCO2_Litev11_Arctic_hist_60E120E.png'])
    saveas(gcf,[path_fig,'SIF_OCO2_Litev11_Arctic_hist_60E120E.fig'])
end

figure % 740 nm dayly
aux_1 = squeeze(D2015_Daily_SIF_740nm(1,4800:6000,:));
aux_2 = squeeze(D2016_Daily_SIF_740nm(1,4800:6000,:));
aux_3 = squeeze(D2017_Daily_SIF_740nm(1,4800:6000,:));
histogram(aux_1(:),'Normalization','pdf','FaceAlpha',0.3,'FaceColor','r')
hold on 
histogram(aux_2(:),'Normalization','pdf','FaceAlpha',0.3,'FaceColor','k')
histogram(aux_3(:),'Normalization','pdf','FaceAlpha',0.3,'FaceColor','g')
legend(['Daily SIF 2015 - mean value:',num2str(round(mean(aux_1(:),'omitnan'),3))],...
       ['Daily SIF 2016 - mean value:',num2str(round(mean(aux_2(:),'omitnan'),3))],...
       ['Daily SIF 2017 - mean value:',num2str(round(mean(aux_3(:),'omitnan'),3))])
set(gca,'Fontsize',18)
grid on 
box on 
xlabel('Daily SIF at 740 nm (mW/(m^2 sr nm))')
ylabel('normalized distribution (pdf)')
title('SIF Area: 60^{\circ}N - 90^{\circ}N, 60^{\circ}E - 120^{\circ}E, ')
xlim([-2,2])
if flag_save
    saveas(gcf,[path_fig,'SIF_Daily740_OCO2_Litev11_Arctic_hist_60E120E.png'])
    saveas(gcf,[path_fig,'SIF_Daily740_OCO2_Litev11_Arctic_hist_60E120E.fig'])
end


figure % 757 nm dayly
aux_1 = squeeze(D2015_Daily_SIF_757nm(1,4800:6000,:));
aux_2 = squeeze(D2016_Daily_SIF_757nm(1,4800:6000,:));
aux_3 = squeeze(D2017_Daily_SIF_757nm(1,4800:6000,:));
histogram(aux_1(:),'Normalization','pdf','FaceAlpha',0.3,'FaceColor','r')
hold on 
histogram(aux_2(:),'Normalization','pdf','FaceAlpha',0.3,'FaceColor','k')
histogram(aux_3(:),'Normalization','pdf','FaceAlpha',0.3,'FaceColor','g')
legend(['Daily SIF 2015 - mean value:',num2str(round(mean(aux_1(:),'omitnan'),3))],...
       ['Daily SIF 2016 - mean value:',num2str(round(mean(aux_2(:),'omitnan'),3))],...
       ['Daily SIF 2017 - mean value:',num2str(round(mean(aux_3(:),'omitnan'),3))])
set(gca,'Fontsize',18)
grid on 
box on 
xlabel('Daily SIF at 757 nm (mW/(m^2 sr nm))')
ylabel('normalized distribution (pdf)')
title('SIF Area: 60^{\circ}N - 90^{\circ}N, 60^{\circ}E - 120^{\circ}E, ')
xlim([-2,2])
if flag_save
    saveas(gcf,[path_fig,'SIF_Daily757_OCO2_Litev11_Arctic_hist_60E120E.png'])
    saveas(gcf,[path_fig,'SIF_Daily757_OCO2_Litev11_Arctic_hist_60E120E.fig'])
end

figure % 771 nm dayly
aux_1 = squeeze(D2015_Daily_SIF_771nm(1,4800:6000,:));
aux_2 = squeeze(D2016_Daily_SIF_771nm(1,4800:6000,:));
aux_3 = squeeze(D2017_Daily_SIF_771nm(1,4800:6000,:));
histogram(aux_1(:),'Normalization','pdf','FaceAlpha',0.3,'FaceColor','r')
hold on 
histogram(aux_2(:),'Normalization','pdf','FaceAlpha',0.3,'FaceColor','k')
histogram(aux_3(:),'Normalization','pdf','FaceAlpha',0.3,'FaceColor','g')
legend(['Daily SIF 2015 - mean value:',num2str(round(mean(aux_1(:),'omitnan'),3))],...
       ['Daily SIF 2016 - mean value:',num2str(round(mean(aux_2(:),'omitnan'),3))],...
       ['Daily SIF 2017 - mean value:',num2str(round(mean(aux_3(:),'omitnan'),3))])
set(gca,'Fontsize',18)
grid on 
box on 
xlabel('Daily SIF at 771 nm (mW/(m^2 sr nm))')
ylabel('normalized distribution (pdf)')
title('SIF Area: 60^{\circ}N - 90^{\circ}N, 60^{\circ}E - 120^{\circ}E')
xlim([-1,1])
if flag_save
    saveas(gcf,[path_fig,'SIF_Daily771_OCO2_Litev11_Arctic_hist_60E120E.png'])
    saveas(gcf,[path_fig,'SIF_Daily771_OCO2_Litev11_Arctic_hist_60E120E.fig'])
end



% ------------------------------------------------------------------
% Histograms Global for the Arctic Region 
% ------------------------------------------------------------------
figure % 740 nm
aux_1 = squeeze(D2015_SIF_740nm(1,:,:));
aux_2 = squeeze(D2016_SIF_740nm(1,:,:));
aux_3 = squeeze(D2017_SIF_740nm(1,:,:));
histogram(aux_1(:),'Normalization','pdf','FaceAlpha',0.3,'FaceColor','r')
hold on 
histogram(aux_2(:),'Normalization','pdf','FaceAlpha',0.3,'FaceColor','k')
histogram(aux_3(:),'Normalization','pdf','FaceAlpha',0.3,'FaceColor','g')
legend(['SIF 2015 - mean value:',num2str(round(mean(aux_1(:),'omitnan'),3))],...
       ['SIF 2016 - mean value:',num2str(round(mean(aux_2(:),'omitnan'),3))],...
       ['SIF 2017 - mean value:',num2str(round(mean(aux_3(:),'omitnan'),3))])
set(gca,'Fontsize',18)
grid on 
box on 
xlabel('SIF at 740 nm (mW/(m^2 sr nm))')
ylabel('normalized distribution (pdf)')
title('SIF Area: 60^{\circ}N - 90^{\circ}N, 180^{\circ}W - 180^{\circ}E')
xlim([-2,2])
if flag_save
    saveas(gcf,[path_fig,'SIF_OCO2_Litev11_Arctic_hist.png'])
    saveas(gcf,[path_fig,'SIF_OCO2_Litev11_Arctic_hist.fig'])
end

figure % 740 nm dayly
aux_1 = squeeze(D2015_Daily_SIF_740nm(1,:,:));
aux_2 = squeeze(D2016_Daily_SIF_740nm(1,:,:));
aux_3 = squeeze(D2017_Daily_SIF_740nm(1,:,:));
histogram(aux_1(:),'Normalization','pdf','FaceAlpha',0.3,'FaceColor','r')
hold on 
histogram(aux_2(:),'Normalization','pdf','FaceAlpha',0.3,'FaceColor','k')
histogram(aux_3(:),'Normalization','pdf','FaceAlpha',0.3,'FaceColor','g')
legend(['Daily SIF 2015 - mean value:',num2str(round(mean(aux_1(:),'omitnan'),3))],...
       ['Daily SIF 2016 - mean value:',num2str(round(mean(aux_2(:),'omitnan'),3))],...
       ['Daily SIF 2017 - mean value:',num2str(round(mean(aux_3(:),'omitnan'),3))])
set(gca,'Fontsize',18)
grid on 
box on 
xlabel('Daily SIF at 740 nm (mW/(m^2 sr nm))')
ylabel('normalized distribution (pdf)')
title('SIF Area: 60^{\circ}N - 90^{\circ}N, 180^{\circ}W - 180^{\circ}E')
xlim([-2,2])
if flag_save
    saveas(gcf,[path_fig,'SIF_Daily740_OCO2_Litev11_Arctic_hist.png'])
    saveas(gcf,[path_fig,'SIF_Daily740_OCO2_Litev11_Arctic_hist.fig'])
end


figure % 757 nm dayly
aux_1 = squeeze(D2015_Daily_SIF_757nm(1,:,:));
aux_2 = squeeze(D2016_Daily_SIF_757nm(1,:,:));
aux_3 = squeeze(D2017_Daily_SIF_757nm(1,:,:));
histogram(aux_1(:),'Normalization','pdf','FaceAlpha',0.3,'FaceColor','r')
hold on 
histogram(aux_2(:),'Normalization','pdf','FaceAlpha',0.3,'FaceColor','k')
histogram(aux_3(:),'Normalization','pdf','FaceAlpha',0.3,'FaceColor','g')
legend(['Daily SIF 2015 - mean value:',num2str(round(mean(aux_1(:),'omitnan'),3))],...
       ['Daily SIF 2016 - mean value:',num2str(round(mean(aux_2(:),'omitnan'),3))],...
       ['Daily SIF 2017 - mean value:',num2str(round(mean(aux_3(:),'omitnan'),3))])
set(gca,'Fontsize',18)
grid on 
box on 
xlabel('Daily SIF at 757 nm (mW/(m^2 sr nm))')
ylabel('normalized distribution (pdf)')
title('SIF Area: 60^{\circ}N - 90^{\circ}N, 180^{\circ}W - 180^{\circ}E')
xlim([-2,2])
if flag_save
    saveas(gcf,[path_fig,'SIF_Daily757_OCO2_Litev11_Arctic_hist.png'])
    saveas(gcf,[path_fig,'SIF_Daily757_OCO2_Litev11_Arctic_hist.fig'])
end

figure % 771 nm dayly
aux_1 = squeeze(D2015_Daily_SIF_771nm(1,:,:));
aux_2 = squeeze(D2016_Daily_SIF_771nm(1,:,:));
aux_3 = squeeze(D2017_Daily_SIF_771nm(1,:,:));
histogram(aux_1(:),'Normalization','pdf','FaceAlpha',0.3,'FaceColor','r')
hold on 
histogram(aux_2(:),'Normalization','pdf','FaceAlpha',0.3,'FaceColor','k')
histogram(aux_3(:),'Normalization','pdf','FaceAlpha',0.3,'FaceColor','g')
legend(['Daily SIF 2015 - mean value:',num2str(round(mean(aux_1(:),'omitnan'),3))],...
       ['Daily SIF 2016 - mean value:',num2str(round(mean(aux_2(:),'omitnan'),3))],...
       ['Daily SIF 2017 - mean value:',num2str(round(mean(aux_3(:),'omitnan'),3))])
set(gca,'Fontsize',18)
grid on 
box on 
xlabel('Daily SIF at 771 nm (mW/(m^2 sr nm))')
ylabel('normalized distribution (pdf)')
title('SIF Area: 60^{\circ}N - 90^{\circ}N, 180^{\circ}W - 180^{\circ}E')
xlim([-1,1])
if flag_save
    saveas(gcf,[path_fig,'SIF_Daily771_OCO2_Litev11_Arctic_hist.png'])
    saveas(gcf,[path_fig,'SIF_Daily771_OCO2_Litev11_Arctic_hist.fig'])
end




% Deltas 
f = figure;
f.Position = [100 100 800 600];
m_proj('stereographic', 'lat',90,'long',30,'radius',30);
m_pcolor(LON,LAT,squeeze(D2016_SIF_740nm(1,:,:)-D2015_SIF_740nm(1,:,:))); % Ensure lon and lat match the grid of land_cover data
shading flat; % Optional: removes grid lines in pcolor plots
m_coast('LineStyle','-','color','k');
m_grid('box','fancy','tickdir','in','fontsize', 18);
h         = colorbar;
% Adjust colorbar position
caxis([-0.3,0.3])
cb_pos    = get(h, 'Position');
cb_pos(1) = cb_pos(1) + 0.05; % Adjust the amount of shift as needed
set(h, 'Position', cb_pos,'Fontsize',18)
ylabel(h,'\Delta SIF @ 740 nm','Fontsize',18)
title('SIF OCO-2 Lite at 740nm v11 2016-2015')
set(gca,'Fontsize',18)



f = figure;
f.Position = [100 100 800 600];
m_proj('stereographic', 'lat',90,'long',30,'radius',30);
m_pcolor(LON,LAT,squeeze(D2017_SIF_740nm(1,:,:)-D2016_SIF_740nm(1,:,:))); % Ensure lon and lat match the grid of land_cover data
shading flat; % Optional: removes grid lines in pcolor plots
m_coast('LineStyle','-','color','k');
m_grid('box','fancy','tickdir','in','fontsize', 18);
h         = colorbar;
% Adjust colorbar position
caxis([-0.3,0.3])
cb_pos    = get(h, 'Position');
cb_pos(1) = cb_pos(1) + 0.05; % Adjust the amount of shift as needed
set(h, 'Position', cb_pos,'Fontsize',18)
ylabel(h,'\Delta SIF @ 740 nm','Fontsize',18)
title('SIF OCO-2 Lite at 740nm v11 2017-2016')
set(gca,'Fontsize',18)