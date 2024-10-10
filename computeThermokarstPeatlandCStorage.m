% ###################################################################
% ################# computeThermokastPeatlandCStorage   ###########
% This script reads the Thermokarst map and the permafrost-affected
% Peatland C storage dataset to check the degree of consistency between
% these two databases. 

% We use the Thermoskast masps available in 
% https://daac.ornl.gov/cgi-bin/dsviewer.pl?ds_id=1332 and
% published in https://www.nature.com/articles/ncomms13043
% to classify the slopes according to the Thermokast presence/degree.

% We use the permafrost-affected peatland C storage from the paper 
% https://www.pnas.org/doi/epdf/10.1073/pnas.1916387117
% to classify the slopes according to the C storage presence/degree.

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
path_in     = config.input_dir.computed;
years_      = 'cSIF_clear';

% add path for the statistic toolbox IPL
addpath(genpath(config.matlab_paths.simpleR_master))

%% ---------------------------
% READ all the datasets 
% ---------------------------

% classification in Thermokast 
% ncdisp(file)
file         = [path_in,'Thermokast_info_all.nc'];
latitude     = ncread(file,'latitude');
longitude    = ncread(file,'longitude');
TSOC_kgC     = ncread(file,'TSOC_kgC'); % Total Soil Organic Carbon in Kg C
% tkwp         = ncread(file,'tkwp'); % Wetland Thermokast terrain coverage
% tkthlp       = ncread(file,'tkthlp'); % Lake Thermokast terrain coverage
% tkhp         = ncread(file,'tkhp'); % Hillslope Thermokast terrain coverage
% tk_all       = ncread(file,'tk_all'); % All Thermokast terrain coverage
CLASS_int    = ncread(file,'tk_all'); % All Thermokast terrain coverage
CLASS_int    = double(CLASS_int)+1; % So the Thermokast all classes are between 1 to 5 to match the j loops


% read the C storage on Peatlands and permafrost-affected peatlands

file         = [path_in,'Histel_SOC_hg_per_sqm_reprojected.nc']; % Histels == permafrost-affected
% ncdisp(file)
% lat     = ncread(file,'lat'); % Already in ArticSIF-study coordinates
% lon     = ncread(file,'lon'); % Already in ArticSIF-study coordinates
Cpeat   = ncread(file,'regridded_data');

% Convert into kgC/m2
Cpeat(Cpeat<0) = NaN;
Cpeat = Cpeat./10;
% Classes none=1, low=2, moderate=3, high=4, and very high=5 
class_none  = Cpeat(CLASS_int==1);
class_low   = Cpeat(CLASS_int==2);
class_mod   = Cpeat(CLASS_int==3);
class_high  = Cpeat(CLASS_int==4);
class_vhigh = Cpeat(CLASS_int==5);


%% ---------------------------
% Ploting figures 
% ---------------------------

figure
subplot(1,5,1)
histogram(class_none(class_none>0),20,'Normalization','probability');
subplot(1,5,2)
histogram(class_low(class_low>0),20,'Normalization','probability');
subplot(1,5,3)
histogram(class_mod(class_mod>0),20,'Normalization','probability');
subplot(1,5,4)
histogram(class_high(class_high>0),20,'Normalization','probability');
subplot(1,5,5)
histogram(class_vhigh(class_vhigh>0),20,'Normalization','probability');


figure
histogram(class_none(class_none>0),200,'Normalization','probability');
hold on
histogram(class_low(class_low>0),200,'Normalization','probability');
histogram(class_mod(class_moc>0),200,'Normalization','probability');
histogram(class_high(class_high>0),200,'Normalization','probability');
histogram(class_vhigh(class_vhigh>0),200,'Normalization','probability');



figure
colors = [227, 227, 227;182, 146, 222;131, 89, 179;110, 46, 184;39, 4, 79]./255;
cmap_edge = 'k';
cmap_face = colors(1,:);
FaceAlpha_val = 0.6;
[ h, AX ] = boxplots(cmap_edge,cmap_face,FaceAlpha_val,class_none(class_none>0),1);
hold on 
cmap_face = colors(2,:);
[ h, AX ] = boxplots(cmap_edge,cmap_face,FaceAlpha_val,class_low(class_low>0),2)
cmap_face = colors(3,:);
[ h, AX ] = boxplots(cmap_edge,cmap_face,FaceAlpha_val,class_mod(class_mod>0),3)
cmap_face = colors(4,:);
[ h, AX ] = boxplots(cmap_edge,cmap_face,FaceAlpha_val,class_high(class_high>0),4)
cmap_face = colors(5,:);
[ h, AX ] = boxplots(cmap_edge,cmap_face,FaceAlpha_val,class_vhigh(class_vhigh>0),5)

set(gca,'Fontsize',18)



% Without outliers 
figure
colors = [227, 227, 227;182, 146, 222;131, 89, 179;110, 46, 184;39, 4, 79]./255;
cmap_edge = colors(1,:);
cmap_face = colors(1,:);
FaceAlpha_val = 0.6;
[ h, AX ] = boxplots(cmap_edge,cmap_face,FaceAlpha_val,class_none,1,'outliers',false);
hold on 
cmap_face = colors(2,:);
cmap_edge = cmap_face;
[ h, AX ] = boxplots(cmap_edge,cmap_face,FaceAlpha_val,class_low,2,'outliers',false);
cmap_face = colors(3,:);
cmap_edge = cmap_face;
[ h, AX ] = boxplots(cmap_edge,cmap_face,FaceAlpha_val,class_mod,3,'outliers',false);
cmap_face = colors(4,:);
cmap_edge = cmap_face;
[ h, AX ] = boxplots(cmap_edge,cmap_face,FaceAlpha_val,class_high,4,'outliers',false);
cmap_face = colors(5,:);
cmap_edge = cmap_face;
[ h, AX ] = boxplots(cmap_edge,cmap_face,FaceAlpha_val,class_vhigh,5,'outliers',false);
set(gca,'Fontsize',18)
tickLabels = {'None','Low','Moderate','High','Very High'}; 
set(gca,'XTickLabel', tickLabels)
box on 
title('Permafrost-affected (Histel) Peatland C storage')
ylabel('Peatland C storage (kg C / m^2)')




