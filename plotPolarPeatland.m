% ###################################################################
% ####### plotPolarPeatland  ####
% This script reads the peatland information (C store stocks in peatlands)
% data and plots the corresponding maps. 
% Original data be found as supplemental data
% for the scientific paper Hugelius et al 2020 (https://www.pnas.org/cgi/doi/10.1073/pnas.1916387117)

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
path_fig = [config.output_dir.general]; % Saving the results in figures
addpath(genpath(config.matlab_paths.mappingToolbox))

% Plotting Saving settings
plotting  = 1;
save_flag = 1;

% Read the NetCDF file 
filename        = [path_in,'Histel_SOC_hg_per_sqm_reprojected.nc']; %Permafrost-affected peatlands
CpeatHistel     = ncread(filename,'regridded_data'); % 'hg / m-2'
lat             = ncread(filename,'lat'); 
lon             = ncread(filename,'lon'); 

filename        = [path_in,'Histosol_SOC_hg_per_sqm_reprojected.nc'];
CpeatHistosol   = ncread(filename,'regridded_data'); % 'hg / m-2'
% Total carbon
Cpeat           = CpeatHistel+CpeatHistosol;


% Read the NDVI NetCDF file to create a land mask
filename = [path_in,'NDVI_EVI_yearly.nc'];
NDVI     = ncread(filename,'NDVI_Arctic');
mask     = NaN([size(NDVI,2),size(NDVI,3)]);
mask(squeeze(NDVI(1,:,:)) > 0) = 1;

% Define color points (RGB values normalized to [0, 1])
color_points = [
    0.8 0.8 0.8;  % Grays
    0 0.5 1;      % Bluish
    0 0.8 0;      % Greenish
    1 1 0;        % Yellowish
    1 0.5 0;      % Orangish
    1 0 0         % Redish
];

% Positions for each color in the colormap (equally spaced from 0 to 1)
positions = linspace(0, 1, size(color_points, 1));

% Number of colors in the colormap
num_colors = 256;  % Adjust as needed for desired colormap resolution

% Interpolate colors to create a smooth colormap
new_map = interp1(positions, color_points, linspace(0, 1, num_colors));


% % V2 maps reprojected by the authors
% % To double-check with the V1 reprojected 
% fil='/Users/neus/Dropbox/IPL/02_PROYECTOS_PROPIOS/SIF_GPP/database_input/Hugelius_etal_2020_PNAS_supplement_grids_10km_Int_WAED/v2/Hugelius_etal_2020_PNAS_grids/NetCDF_WGS84/Histel_SOC_hg_per_sqm_WGS84.nc';
% ncdisp(fil)
% aa = ncread(fil,'Histel_SOC_hg_per_sqm_WGS84.tif');
% lat_aa = ncread(fil,'lat');
% lon_aa = ncread(fil,'lon');
% % Plotting Estimated storage of organic carbon in peat of Histel peatlands (unit is hg / m-2)'
% f = figure;
% f.Position = [100 100 900 600];
% m_proj('stereographic', 'lat',90,'long',30,'radius',30);
% aux = squeeze(CpeatHistel);
% aux(isnan(mask))=NaN; % Non valid
% m_pcolor(lon_aa,lat_aa,aa');
% shading flat; % Optional: removes grid lines in pcolor plots
% m_coast('LineStyle','-','color','k');
% m_grid('box','fancy','tickdir','in','fontsize', 18);
% h = colorbar;
% % Adjust colorbar position
% cb_pos    = get(h, 'Position');
% cb_pos(1) = cb_pos(1) + 0.05; % Adjust the amount of shift as needed
% set(h, 'Position', cb_pos,'Fontsize',18)
% colormap(new_map)
% ylabel(h,'hg C m^{-2}','Fontsize',18)
% t = title(['Permafrost-affected peatland C storage Author']);
% set(gca,'Fontsize',18)
% set(t,'Position',get(t,'Position')+[0 .04 0]);  % move up slightly
% caxis([0,300])



if plotting == 1
    % Plotting Estimated storage of organic carbon in peat of Histel peatlands (unit is hg / m-2)'
    f = figure;
    f.Position = [100 100 900 600];    
    m_proj('stereographic', 'lat',90,'long',30,'radius',30);
    aux = squeeze(CpeatHistel);
    aux(isnan(mask))=NaN; % Non valid
    m_pcolor(lon,lat,aux');
    shading flat; % Optional: removes grid lines in pcolor plots
    m_coast('LineStyle','-','color','k');
    m_grid('box','fancy','tickdir','in','fontsize', 22);
    h = colorbar;
    % Adjust colorbar position
    cb_pos    = get(h, 'Position');
    cb_pos(1) = cb_pos(1) + 0.05; % Adjust the amount of shift as needed
    set(h, 'Position', cb_pos,'Fontsize',22)
    colormap(new_map)
    ticklabels = {'0', '50', '100', '150', '200','>300'}; % Define colorbar tick labels
    caxis([0, 300]);  % Extend the colorbar range to include values up to 350
    ticks = linspace(0, 300, length(ticklabels));
    set(h, 'Ticks', ticks, 'TickLabels', ticklabels);
    ylabel(h, 'hg C m^{-2}', 'FontSize', 22);
    t = title(['Permafrost-affected peatland C storage']);
    set(gca,'Fontsize',18)
    set(t,'Position',get(t,'Position')+[0 .04 0]);  % move up slightly

    
    if save_flag
        saveas(gcf,[path_fig,'PermafrostPeatlandC_025deg.png'])
        saveas(gcf,[path_fig,'PermafrostPeatlandC_025deg.fig'])
    end
    
     % Plotting Estimated storage of organic carbon in peat of Histel peatlands (unit is hg / m-2)'
    f = figure;
    f.Position = [100 100 900 600];    
    m_proj('stereographic', 'lat',90,'long',30,'radius',30);
    aux = squeeze(CpeatHistosol);
    aux(isnan(mask))=NaN; % Non valid
    m_pcolor(lon,lat,aux');
    shading flat; % Optional: removes grid lines in pcolor plots
    m_coast('LineStyle','-','color','k');
    m_grid('box','fancy','tickdir','in','fontsize', 18);
    h = colorbar;
    % Adjust colorbar position
    cb_pos    = get(h, 'Position');
    cb_pos(1) = cb_pos(1) + 0.05; % Adjust the amount of shift as needed
    set(h, 'Position', cb_pos,'Fontsize',18)
    colormap(new_map)
    ticklabels = {'0', '50', '100', '150', '200','>300'}; % Define colorbar tick labels
    caxis([0, 300]);  % Extend the colorbar range to include values up to 350
    ticks = linspace(0, 300, length(ticklabels));
    set(h, 'Ticks', ticks, 'TickLabels', ticklabels);
    ylabel(h,'hg C m^{-2}','Fontsize',18)
    t = title(['Peatland C storage']);
    set(gca,'Fontsize',18)
    set(t,'Position',get(t,'Position')+[0 .04 0]);  % move up slightly
    
    if save_flag
        saveas(gcf,[path_fig,'PeatlandC_025deg.png'])
        saveas(gcf,[path_fig,'PeatlandC_025deg.fig'])
    end
    
     % Plotting Estimated storage of organic carbon in peat of Histel peatlands (unit is hg / m-2)'
    f = figure;
    f.Position = [100 100 900 600];    
    m_proj('stereographic', 'lat',90,'long',30,'radius',30);
    aux = squeeze(Cpeat);
    aux(isnan(mask))=NaN; % Non valid
    m_pcolor(lon,lat,aux');
    shading flat; % Optional: removes grid lines in pcolor plots
    m_coast('LineStyle','-','color','k');
    m_grid('box','fancy','tickdir','in','fontsize', 18);
    h = colorbar;
    % Adjust colorbar position
    cb_pos    = get(h, 'Position');
    cb_pos(1) = cb_pos(1) + 0.05; % Adjust the amount of shift as needed
    set(h, 'Position', cb_pos,'Fontsize',18)
    colormap(new_map)
    ticklabels = {'0', '50', '100', '150', '200','>300'}; % Define colorbar tick labels
    caxis([0, 300]);  % Extend the colorbar range to include values up to 350
    ticks = linspace(0, 300, length(ticklabels));
    set(h, 'Ticks', ticks, 'TickLabels', ticklabels);
    ylabel(h,'hg C m^{-2}','Fontsize',18)
    t = title(['Total Peatland C storage']);
    set(gca,'Fontsize',18)
    set(t,'Position',get(t,'Position')+[0 .04 0]);  % move up slightly
    
    if save_flag
        saveas(gcf,[path_fig,'TotalPeatlandC_025deg.png'])
        saveas(gcf,[path_fig,'TotalPeatlandC_025deg.fig'])
    end
end
