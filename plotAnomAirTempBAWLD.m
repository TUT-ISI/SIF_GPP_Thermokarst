% ###################################################################
% ##################### PlotAnoAirTempBAWLD  ########################
% ###################################################################

% Plot Air Temperature anomalies saved and estimated in the period of
% time of interest for the study and evaluate its range per class of BAWLD


% Author: Neus Sabater
% e-mail: neus.sabater@fmi.fi 
% version 0: (19/03/2024) 
% --------------------------

% Main script
addpath(genpath([pwd,filesep,'functions_environment']))
% Automatically detect environment
env = detect_environment();
% Load configuration
config = load_conf(env);

%% --------------------------
%  set parameters 
% ---------------------------

plotting    = 1;
save_flag   = 0;

path_fig    = config.output_dir.general; %'/Users/neus/Dropbox/IPL/02_PROYECTOS_PROPIOS/SIF_GPP/figures_all/figs/';
path_in     = config.inpinput_dir.computed; %'/Users/neus/Dropbox/IPL/02_PROYECTOS_PROPIOS/SIF_GPP/codes/data/Computed_data/';
years_      = 'cSIF_clear';

% Air Temperature anomaly 
load([path_in,'airTAnomaly_JunJul_005deg_2003_2018.mat'])


% Classification in BAWLD
CLASS_int = ncread([path_in,'BAWLD_dataset.nc'],'CLASS_int');

LAT = ncread([path_in,'BAWLD_dataset.nc'],'LAT');
LON = ncread([path_in,'BAWLD_dataset.nc'],'LON');


% Store the temperature per class
n = 15; % Number of classes
for j = 1:n
    aux = []; AUX = [];
    for i = 1:size(t2m_anomaly_mean,3) % years
        t2m = squeeze(t2m_anomaly_mean(:,:,i));
        aux = t2m(CLASS_int==j);
        aux = aux(:);
        AUX = [AUX;aux];        
    end
    eval(['t2m_',num2str(j),'class = AUX;'])
end




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



figure
x = [-5:0.1:5];
for j=1:n
    subplot(4,4,j)
    eval(['H = histogram(t2m_',num2str(j),'class,''Normalization'',''probability'',''BinWidth'',0.1,''FaceColor'',colors(j,:))'])
    mean_vals(j) = mean(H.Data,'omitnan');
    std_vals(j)  = std(H.Data,'omitnan');
%     y     = 1/(sqrt(2*pi)*std_vals(j)).*exp(-(x-mean_vals(j)).^2)./(2.*std_vals(j).^2);
%     hold on 
%     plot(x,y,'color',colors(j,:))
    hold on 
    title(tickLabels{j})
    set(gca,'Fontsize',18)
    xlim([-5,5])
    box on 
    grid on 
    xlabel('Air Temperature anomaly [\circ C]')
    ylabel('Probability [-]')
end 

if save_flag
    saveas(gcf,[path_fig,'BAWLD_airTAnom_hist_subplot.png'])
    saveas(gcf,[path_fig,'BAWLD_airTAnom_hist_subplot.fig'])
end

figure

for j=1:n
    param = [std_vals(j) mean_vals(j)];
    y     = exp(-(x-mean_vals(j)).^2)./(2.*std_vals(j).^2);
    plot(x,y,'color',colors(j,:))
    hold on 
end 


figure
for j=1:n
    eval(['histogram(t2m_',num2str(j),'class,''Normalization'',''probability'',''BinWidth'',0.1,''FaceColor'',colors(j,:),''FaceAlpha'',0.3)'])
    set(gca,'Fontsize',18)
    hold on 
    xlim([-5,5])
end
box on
grid on
xlabel('Air Temperature anomaly [\circ C]')
ylabel('Probability [-]')
legend(tickLabels)
title('Normalised probability of air temperature anomalies')
if save_flag
    saveas(gcf,[path_fig,'BAWLD_airTAnom_hist.png'])
    saveas(gcf,[path_fig,'BAWLD_airTAnom_hist.fig'])
end

