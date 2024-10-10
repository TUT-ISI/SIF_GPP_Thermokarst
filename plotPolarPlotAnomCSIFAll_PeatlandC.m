% ###################################################################
% ################# plotPolarPlotAnomCSIFAll_Thermokast   ###########
% This script reads the cSIF database anomalies from 2003-2016 and computes 
% the OLS regression with the corresponding airT anomalies for the same
% years.  It can also read the cSIF database for the clear computations
% from 2003-2018
% We also use the permafrost-affected peatland C storage from the paper 
% https://www.pnas.org/doi/epdf/10.1073/pnas.1916387117
% to classify the slopes according to the C storage presence/degree.

% Data was extracted from the cSIF database generated with OCO-2 SIF data and 
% MODIS surface reflectance every 4-day at 0.05 degrees. 
% This database also disentangles the sunlit and all sky SIF prediction.
% 
% %Zhang, Y., Joiner, J., Alemohammad, S. H., Zhou, S., & Gentine, P. (2018).
% A global spatially contiguous solar-induced fluorescence (CSIF) dataset 
% using neural networks. Biogeosciences, 15(19), 5779-5800.
% https://bg.copernicus.org/articles/15/5779/2018/

% Author: Neus Sabater
% e-mail: neus.sabater@fmi.fi 
% version 0: (27/06/2024) 
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
path_in     = config.input_dir.computed;
years_      = 'cSIF_clear';

% add path for the statistic toolbox IPL
addpath(genpath(config.matlab_paths.simpleR_master))



%% ---------------------------
% READ all the datasets 
% ---------------------------
% airTemperatureswitch years_ % Changed for the ARCLIM Temperature anomalies
switch years_
    case 'cSIF_all'
%         load([path_in,'airTAnomaly_JunJul_005deg_2003_2016.mat']) % 'years_cov','t2m_summer_anomaly_year','lon','lat'))
%         t2m_anomaly_mean = permute(t2m_anomaly_mean,[2,3,1]);     % change the dimensions long lat years
        
    case 'cSIF_clear'
%         load([path_in,'airTAnomaly_JunJul_005deg_2003_2018.mat']) % 'years_cov','t2m_summer_anomaly_year','lon','lat')) 
        %ncdisp([path_in,'ERA5_ARCLIM_mean_anomalies_2003_2018.nc'])
        t2m_anomaly_mean = ncread([path_in,'ERA5_ARCLIM_mean_anomalies_2003_2018.nc'],'t2m_anom_r');
        latitude         = ncread([path_in,'ERA5_ARCLIM_mean_anomalies_2003_2018.nc'],'latitude');
        longitude        = ncread([path_in,'ERA5_ARCLIM_mean_anomalies_2003_2018.nc'],'longitude');
        t2m_mean         = ncread([path_in,'ERA5_ARCLIM_mean_anomalies_2003_2018.nc'],'t2m__mean_r');
        t2m_anomaly_mean = permute(t2m_anomaly_mean,[2,3,1]); 
end

% read the C storage on Peatlands and permafrost-affected peatlands

file         = [path_in,'Histel_SOC_hg_per_sqm_reprojected.nc']; % permafrost-affected
% ncdisp(file)
% lat     = ncread(file,'lat'); % Already in ArticSIF-study coordinates
% lon     = ncread(file,'lon'); % Already in ArticSIF-study coordinates
Cpeat   = ncread(file,'regridded_data');



% GPP 
load([path_in,'SIF_GPP_Z_score_TrendArticSIFGPPFluxSat_JJ_2003_2018.mat'],'GPP_TrendStore')  
switch years_
    case 'cSIF_all'
        GPP_TrendStore     = GPP_TrendStore(:,:,1:14);% trim the dataseries up to 2016
end
GPP_mean           = mean(GPP_TrendStore,3,'omitnan');
GPP_median         = median(GPP_TrendStore,3,'omitnan');
GPP_anomaly_mean   = GPP_TrendStore-GPP_mean;
GPP_anomaly_median = GPP_TrendStore-GPP_median;



switch years_
    
    case 'cSIF_all'
        % cSIF data
        load([path_in,'cSIFAnomaly_JunJul_05deg_2003_2016.mat'],'cSIF_anomaly_mean','cSIF_anomaly_median',...
            'cSIF_median','cSIF_mean','longitude','latitude','years_covered_cSIF')
        cSIF_anomaly_mean  = permute(cSIF_anomaly_mean,[2,3,1]);   % change the dimensions long lat years
        
        
    case 'cSIF_clear'
        load([path_in,'cSIFclearAnomaly_JunJul_005deg_2003_2018.mat'],'cSIF_anomaly_mean','cSIF_anomaly_median',...
            'cSIF_median','cSIF_mean','longitude','latitude','years_covered_cSIF')
        cSIF_anomaly_mean  = permute(cSIF_anomaly_mean,[2,3,1]);   % change the dimensions long lat years
end

%% --------------------------
% Plotting and computing the slope
% ---------------------------


% if plotting == 1
% GPP anomalies scatter plots vs.C content on Peatlands
f = figure;
f.Position = [100 100 2000 1200];

auxGPP_an_store = [];
auxT_an_store   = [];

for i=1:size(GPP_TrendStore,3) % each year
    auxGPP = squeeze(GPP_anomaly_mean(:,:,i));
    auxT   = squeeze(t2m_anomaly_mean(:,:,i));        
    auxGPP_an_store = [auxGPP_an_store(:);auxGPP(:)]; % all years
    auxT_an_store   = [auxT_an_store(:);auxT(:)]; % all years
end
auxPeatland_an_store = repmat(Cpeat(:),size(GPP_TrendStore,3),1);


% d = scatter(auxPeatland_an_store(auxPeatland_an_store>0),auxGPP_an_store(auxPeatland_an_store>0),8,'filled');
% scatter3(auxPeatland_an_store(auxPeatland_an_store>0)',auxGPP_an_store(auxPeatland_an_store>0)',auxT_an_store(auxPeatland_an_store>0)');


% Example vectors
x = double(auxPeatland_an_store(auxPeatland_an_store>0));
y = auxGPP_an_store(auxPeatland_an_store>0);
z = auxT_an_store(auxPeatland_an_store>0);



% Create a grid
[xq, yq] = meshgrid(linspace(min(x), max(x), 100), linspace(min(y), max(y), 100));
xq = double(xq);
yq = double(yq);
% remove NaN
x = x(~isnan(y) & ~isnan(z));
y_new = y(~isnan(y) & ~isnan(z));
z_new =double(z(~isnan(y) & ~isnan(z)));


% Interpolate the scattered data onto the grid
zq = griddata(x, y_new, z_new, xq, yq, 'linear');

% % Plot the contour
% contour3(zq,yq,xq,20); % 20 is the number of contour levels
% zlabel('Permaforst-affetced Peatland C storage (kgC/m^2)');
% ylabel('GPP anomaly (kgC/m^2)');
% xlabel('Aux T anomaly');
% grid on;
% title('3D Contour Plot');

%% 
figure
s = mesh(xq./10,yq,zq)
 s.FaceColor = 'inter';
set(gca,'Fontsize',18)
xlabel('Permaforst-affetced Peatland C storage (kgC/m^2)')
ylabel('GPP anomaly (kgC/m^2)')
view(0, 90);
box on 
h = colorbar;
caxis([-4,4])
ylabel(h,'air T anomaly (C^{\circ})')
cmap = redblue;
colormap(cmap)
box on
grid on
grid minor
set(gca,'TickLength',[0.05, 0.01])
set(gca,'XMinorTick','on','YMinorTick','on')



if save_flag
    switch years_
        case 'initial'
            saveas(gcf,[path_fig,'PeatlandC_GPP_class_scatter_2002_2016.png'])
            saveas(gcf,[path_fig,'PeatlandC_GPP_class_scatter_2002_2016.fig'])
            
        case 'final'
            saveas(gcf,[path_fig,'PeatlandC_GPP_class_scatter_2003_2018.png'])
            saveas(gcf,[path_fig,'PeatlandC_GPP_class_scatter_2003_2018.fig'])
            
        case 'cSIF_all'
            saveas(gcf,[path_fig,'PeatlandC_GPP_class_scatter_2003_2016.png'])
            saveas(gcf,[path_fig,'PeatlandC_GPP_class_scatter_2003_2016.fig'])
            
        case 'cSIF_clear'
            saveas(gcf,[path_fig,'PeatlandC_GPPclear_class_scatter_2003_2018.png'])
            saveas(gcf,[path_fig,'PeatlandC_GPPclear_class_scatter_2003_2018.fig'])
            
    end
end

    
    

