% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% PlotPolarThermokast: This script read the maps of thermokast of permafrost 
% available in https://daac.ornl.gov/cgi-bin/dsviewer.pl?ds_id=1332 and
% published in https://www.nature.com/articles/ncomms13043
% The original maps (in *shp and *prj formats) have been regridded to the
% lat lon grid used in the ArcticSIF study. The function to generate the
% regridded maps and store them in NetCDF is readThermokast.py

% Version: 0
% data: Mar/2024
% author: neus sabater
% e-mail: neus.sabater@fmi.fi
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


addpath(genpath([pwd,filesep,'functions_environment']))
% Automatically detect environment
env = detect_environment();
% Load configuration
config = load_conf(env);

path_in      = config.input_dir.computed; 
file         = [path_in,'Thermokast_info_all.nc'];
path_fig     = config.output_dir.genearl; 
save_flag    = 0; % 1 active


% ncdisp(file)
latitude     = ncread(file,'latitude');
longitude    = ncread(file,'longitude');
TSOC_kgC     = ncread(file,'TSOC_kgC'); % Total Soil Organic Carbon in Kg C
tkwp         = ncread(file,'tkwp'); % Wetland thermokarst terrain coverage
tkthlp       = ncread(file,'tkthlp'); % Lake thermokarst terrain coverage
tkhp         = ncread(file,'tkhp'); % Hillslope thermokarst terrain coverage
tk_all       = ncread(file,'tk_all'); % All thermokarst terrain coverage


% Example of Total Soil Organic Carbon in Kg C
f = figure;
f.Position = [100 100 700 600];          
m_proj('stereographic', 'lat',90,'long',30,'radius',30);
m_pcolor(double(longitude), double(latitude),double(TSOC_kgC)')
shading flat;
c = colorbar;
% Adjust colorbar position
set(c, 'Location', 'southoutside','Fontsize',18)
ylabel(c,'Total Soil Organic Carbon in Kg C')
m_coast('color', 'k');
m_grid('box', 'fancy', 'fontsize', 16);
t = title('Total Soil Organic Carbon in Kg C');
set(t,'Position',get(t,'Position')+[0 .03 0],'Fontsize',18);  % move up slightly
set(gca,'Fontsize',18)
caxis([0,7e10])

if save_flag
    saveas(gcf,[path_fig,'TotalOrgCarSoil.png'])
    saveas(gcf,[path_fig,'TotalOrgCarSoil.fig'])
end 



% Example of Wetland thermokarst terrain coverage
f = figure;
f.Position = [100 100 700 600];          
m_proj('stereographic', 'lat',90,'long',30,'radius',30);
m_pcolor(double(longitude), double(latitude),double(tkwp)')
shading flat;
c = colorbar;
% Adjust colorbar position
set(c, 'Location', 'southoutside','Fontsize',18)
cmap = [227, 227, 227;177, 240, 182;118, 232, 126;70, 179, 77;25, 74, 28]./255;
colormap(cmap);
m_coast('color', 'k');
m_grid('box', 'fancy', 'fontsize', 16);
c.Ticks = [0.4, 1.2, 2, 2.8, 3.6]; % Positions of ticks
c.TickLabels = {'None','Low','Moderate','High','Very High'}; 
t = title('Wetland thermokarst terrain coverage');
set(t,'Position',get(t,'Position')+[0 .03 0],'Fontsize',18);  % move up slightly

if save_flag
    saveas(gcf,[path_fig,'WetThermokast.png'])
    saveas(gcf,[path_fig,'WetThermokast.fig'])
end 


% Example of Lake thermokarst terrain coverage
f = figure;
f.Position = [100 100 700 600];          
m_proj('stereographic', 'lat',90,'long',30,'radius',30);
m_pcolor(double(longitude), double(latitude),double(tkthlp)')
shading flat;
c = colorbar;
% Adjust colorbar position
set(c, 'Location', 'southoutside','Fontsize',18)
cmap = [227, 227, 227;75, 206, 242;53, 138, 242;28, 86, 158;19, 41, 69]./255;
colormap(cmap);
m_coast('color', 'k');
m_grid('box', 'fancy', 'fontsize', 16);
caxis([0,4])
c.Ticks = [0.4, 1.2, 2, 2.8, 3.6]; % Positions of ticks
c.TickLabels = {'None','Low','Moderate','High','Very High'}; 
t = title('Lake thermokarst terrain coverage');
set(t,'Position',get(t,'Position')+[0 .03 0],'Fontsize',18);  % move up slightly

if save_flag
    saveas(gcf,[path_fig,'LakeThermokast.png'])
    saveas(gcf,[path_fig,'LakeThermokast.fig'])
end 


% Example of Hill thermokarst terrain coverage
f = figure;
f.Position = [100 100 700 600];          
m_proj('stereographic', 'lat',90,'long',30,'radius',30);
m_pcolor(double(longitude), double(latitude),double(tkhp)')
shading flat;
c = colorbar;
% Adjust colorbar position
set(c, 'Location', 'southoutside','Fontsize',18)
cmap = [227, 227, 227;242, 160, 160;237, 104, 104;186, 50, 50;120, 12, 12]./255;
colormap(cmap);
m_coast('color', 'k');
m_grid('box', 'fancy', 'fontsize', 16);
caxis([0,4])
c.Ticks = [0.4, 1.2, 2, 2.8, 3.6]; % Positions of ticks
c.TickLabels = {'None','Low','Moderate','High','Very High'}; 
t = title('Hillslope terrain thermokarst terrain coverage');
set(t,'Position',get(t,'Position')+[0 .03 0],'Fontsize',18);  % move up slightly

if save_flag
    saveas(gcf,[path_fig,'HillThermokast.png'])
    saveas(gcf,[path_fig,'HillThermokast.fig'])
end 


% Example of thermokarst ALL coverage
f = figure;
f.Position = [100 100 700 600];          
m_proj('stereographic', 'lat',90,'long',30,'radius',30);
m_pcolor(double(longitude), double(latitude),double(tk_all)')
shading flat;
c = colorbar;
% Adjust colorbar position
set(c, 'Location', 'southoutside','Fontsize',18)
cmap = [227, 227, 227;182, 146, 222;131, 89, 179;110, 46, 184;39, 4, 79]./255;
colormap(cmap);
m_coast('color', 'k');
m_grid('box', 'fancy', 'fontsize', 16);
c.Ticks = [0.4, 1.2, 2, 2.8, 3.6]; % Positions of ticks
c.TickLabels = {'None','Low','Moderate','High','Very High'}; 
t = title('All thermokarst terrain coverage');
set(t,'Position',get(t,'Position')+[0 .03 0],'Fontsize',18);  % move up slightly

if save_flag
    saveas(gcf,[path_fig,'AllThermokast.png'])
    saveas(gcf,[path_fig,'AllThermokast.fig'])
end 



