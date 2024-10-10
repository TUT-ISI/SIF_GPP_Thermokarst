% ###################################################################
% ################# plotConfusionMatrixThermokarstBAWLD   ###########
% This script reads the Thermokarst classification and the BAWLD class
% and generates a kind matrix that relates the percent of the Thermokarst
% coverage for each class.
% Thermokarst data:
% https://daac.ornl.gov/cgi-bin/dsviewer.pl?ds_id=1332 and
% published in https://www.nature.com/articles/ncomms13043
% to classify the slopes according to the Thermokast presence/degree.

% BAWLD data:
% lake arctic database published here
% https://arcticdata.io/catalog/view/doi:10.18739/A2C824F9X

% Author: Neus Sabater
% e-mail: neus.sabater@fmi.fi 
% version 0: (19/03/2024) 
%% ---------------------------
%  set parameters 
% ----------------------------
plotting    = 1;
save_flag   = 0;
addpath(genpath([pwd,filesep,'functions_environment']))
% Automatically detect environment
env = detect_environment();
% Load configuration
config = load_conf(env);

path_fig    = config.output_dir.general; 
path_in     = config.input_dir.computed; 


% Thermokarst
file         = [path_in,'Thermokast_info_all.nc'];
latitude     = ncread(file,'latitude');
longitude    = ncread(file,'longitude');
CLASS_th     = ncread(file,'tk_all'); % All Thermokast terrain coverage
CLASS_th     = double(CLASS_th)+1; % So the Thermokast all classes are between 1 to 5 to match the j loops
n_th         = 5; % Number of colors according to the Thermokast classification
% Define custom colormap with high contrast
colors_th    = [227, 227, 227;182, 146, 222;131, 89, 179;110, 46, 184;39, 4, 79]./255;
tickLabels_th= {'None','Low','Moderate','High','Very High'};




% Classification in BAWLD
CLASS_int = ncread([path_in,'BAWLD_dataset.nc'],'CLASS_int');
LAT       = ncread([path_in,'BAWLD_dataset.nc'],'LAT');
LON       = ncread([path_in,'BAWLD_dataset.nc'],'LON');
% There are 15 different classes
n_cl      = 15; % Number of colors
% Define custom colormap with high contrast
colors_cl = [
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
colors_cl = colors_cl(1:n_cl, :);
colors_cl = colors_cl ./ max(colors_cl, [], 1);

tickLabels_cl = {'Permafrost Peatlands','Sparse Boreal Peatlands','Rivers',...
    'Glaciers','Upland Tundra','Common Boreal Peatlands','Large Lakes',...
    'Lake-rich wetlands','Dominant Boreal Peatlands','Wetland-rich Tundra',...
    'Alpine and Tundra Barrens','Wetland and Lake-rich Tundra','Lake-rich Shield','Upland Boreal','Wetland and Lake-rich Yedoma Tundra'};



% Create a matrix 
% x-axis the classes y-axis the thermokarst coverage
DATA  = NaN(n_cl,n_th);
count = ones(size(CLASS_int));
for i = 1:n_cl
    num_classi = sum(count(CLASS_int==i)); % Total pixels of class i
    for j = 1:n_th 
        
        num_classj = sum(count(CLASS_int==i & CLASS_th==j));
        DATA(i,j)  = 100.*(num_classj./num_classi);
        
    end 
end 




data = DATA; 
figure
for i = 1:size(data, 1)
    for j = 1:size(data, 2)
        alpha_val = data(i, j)./100; % Use the data percent value as alpha
        alpha_val = 0.2+alpha_val; % Sum 0.2 for visibility
        h = rectangle('Position', [j-0.5, i-0.5, 1, 1], 'FaceColor', [colors_cl(i, :), alpha_val], 'EdgeColor', 'none');
    end
end

box on 
ylim([0.5,15.5])
xlim([0.5,5.5])
xticks([1:5])
xticklabels({'None','Low','Moderate','High','Very High'}); 
yticks([1:15])
yticklabels({'Permafrost Peatlands','Sparse Boreal Peatlands','Rivers',...
    'Glaciers','Upland Tundra','Common Boreal Peatlands','Large Lakes',...
    'Lake-rich wetlands','Dominant Boreal Peatlands','Wetland-rich Tundra',...
    'Alpine and Tundra Barrens','Wetland and Lake-rich Tundra','Lake-rich Shield','Upland Boreal','Wetland and Lake-rich Yedoma Tundra'});

set(gca,'Fontsize',18)
xlabel('Thermokarst coverage (%)')
ylabel('BAWLD classification')
% Create a rectangle around the entire matrix
box_position = [0.5, 0.5, size(data, 2), size(data, 1)];
rectangle('Position', box_position, 'EdgeColor', 'k', 'LineWidth', 2);
% Remove x-axis and y-axis ticks
set(gca, 'TickLength', [0, 0]); % Set y-axis tick length to 0
hold on 
for i = 1:size(data, 1)  
    plot((0.5:15.5),(i+0.5).*ones(1,16),'k-')
end


% Customize colorbar labels
c = colorbar;
colormap(gray)
c.Label.String = 'Percent of Thermokarst (%)';
c.Ticks = [0, 1];
c.TickLabels = {'High', 'Low'};
c.FontSize = 18;


if save_flag
    saveas(gcf,[path_fig,'Thermokast_BAWLD_percentmatrix.png'])
    saveas(gcf,[path_fig,'Thermokast_BAWLD_percentmatrix.fig'])
end


