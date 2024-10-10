% ###################################################################
% ####### Percents in  BAWLD ####
% This script reads the BAWLD database together with the percent
% information to compute a stacked bar chart to estimate the percent of
% different elements within the grid cells classified as each class.
% Percents of:
% Percent of galciers and permanent snow within grid cell (%).  
% Percent of rockland within grid cell (%).
% Percent of tundra within grid cell (%).
% Percent of boreal forest within grid cell (%).
% Percent of wetland within grid cell (%) Includes, permafrost bogs, wet tundra, marsh, bog, and fen.
% Percent of lakes within grid cell (%).
% Percent of rivers within grid cell (%).
% K-means classification. 

% lake arctic database published here
% https://arcticdata.io/catalog/view/doi:10.18739/A2C824F9X

% Reading routine of the *.shp files and NetCDF storage in readBAWLD.m

% Author: Neus Sabater
% e-mail: neus.sabater@fmi.fi 
% version 0: (19/03/2024) 
%% ---------------------------
%  set parameters 
% ---------------------------
addpath(genpath([pwd,filesep,'functions_environment']))

% Automatically detect environment
env = detect_environment();

% Load configuration
config = load_conf(env);

plotting    = 1;
save_flag   = 0;

path_in     = config.input_dir.computed; 
path_fig    = config.output_dir.general; 
years_      = 'cSIF_clear';



% Classification in BAWLD
CLASS_int = ncread([path_in,'BAWLD_dataset.nc'],'CLASS_int');
% Components 
BOR_int = ncread([path_in,'BAWLD_dataset.nc'],'BOR_int');
GLA_int = ncread([path_in,'BAWLD_dataset.nc'],'GLA_int');
TUN_int = ncread([path_in,'BAWLD_dataset.nc'],'TUN_int');
WET_int = ncread([path_in,'BAWLD_dataset.nc'],'WET_int');
LAK_int = ncread([path_in,'BAWLD_dataset.nc'],'LAK_int');
RIV_int = ncread([path_in,'BAWLD_dataset.nc'],'RIV_int');
ROC_int = ncread([path_in,'BAWLD_dataset.nc'],'ROC_int');
LAT = ncread([path_in,'BAWLD_dataset.nc'],'LAT');
LON = ncread([path_in,'BAWLD_dataset.nc'],'LON');


% There are 15 different classes
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

tickLabels = {'Permafrost Peatlands','Sparse Boreal Peatlands','Rivers',...
    'Glaciers','Upland Tundra','Common Boreal Peatlands','Large Lakes',...
    'Lake-rich wetlands','Dominant Boreal Peatlands','Wetland-rich Tundra',...
    'Alpine and Tundra Barrens','Wetland and Lake-rich Tundra','Lake-rich Shield','Upland Boreal','Wetland and Lake-rich Yedoma Tundra'};


for i =1:n
    BOR_aux = BOR_int(CLASS_int==i);
    BOR_mean = mean(BOR_aux(:),1,'omitnan');
    
    GLA_aux = GLA_int(CLASS_int==i);
    GLA_mean = mean(GLA_int(:),1,'omitnan');
    
    TUN_aux = TUN_int(CLASS_int==i);
    TUN_mean = mean(TUN_aux(:),1,'omitnan');
    
    WET_aux = WET_int(CLASS_int==i);
    WET_mean = mean(WET_aux(:),1,'omitnan');
    
    LAK_aux = LAK_int(CLASS_int==i);
    LAK_mean = mean(LAK_aux(:),1,'omitnan');
    
    RIV_aux = RIV_int(CLASS_int==i);
    RIV_mean = mean(RIV_aux(:),1,'omitnan');
    
    ROC_aux = ROC_int(CLASS_int==i);
    ROC_mean = mean(ROC_aux(:),1,'omitnan');
    
    percent_class(i,:) = [BOR_mean,GLA_mean,TUN_mean,WET_mean,LAK_mean,RIV_mean,ROC_mean];
        
end 


figure
barh(percent_class,'stacked')
legend('boreal forest','galciers','tundra','wetland','lakes','rivers','rockland')
set(gca,'Fontsize',18)
yticklabels({'Permafrost Peatlands','Sparse Boreal Peatlands','Rivers',...
    'Glaciers','Upland Tundra','Common Boreal Peatlands','Large Lakes',...
    'Lake-rich wetlands','Dominant Boreal Peatlands','Wetland-rich Tundra',...
    'Alpine and Tundra Barrens','Wetland and Lake-rich Tundra','Lake-rich Shield','Upland Boreal','Wetland and Lake-rich Yedoma Tundra'});
% xtickangle(45)
set(get(gca, 'XAxis'), 'FontWeight', 'bold');
box on
grid on
grid minor
xlim([0,110])
xlabel('Percent (%)')

title('BAWLD database')

saveas(gcf,[path_fig,'BAWLD_percent_horizontal.png'])
saveas(gcf,[path_fig,'BAWLD_percent_horizontal.fig'])


