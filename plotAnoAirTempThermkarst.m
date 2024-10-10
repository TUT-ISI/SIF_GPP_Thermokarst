% ###################################################################
% ##################### PlotAnoAirTempThemokarst  ###################
% ###################################################################

% Plot air Temperature anomalies saved and estimated in the period of
% time of interest for the study and evaluate its range per class of BAWLD

% Author: Neus Sabater
% e-mail: neus.sabater@fmi.fi 
% version 0: (19/03/2024) 
% ---------------------------

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
path_in     = config.input_dir.computed; %'/Users/neus/Dropbox/IPL/02_PROYECTOS_PROPIOS/SIF_GPP/codes/data/Computed_data/';
years_      = 'cSIF_clear';

% Air Temperature anomaly 
load([path_in,'airTAnomaly_JunJul_005deg_2003_2018.mat'])


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

n_years = size(t2m_anomaly_mean,3);
AUX_1 = [];
AUX_2 = [];
AUX_3 = [];
AUX_4 = [];
AUX_5 = [];

for i = 1:n_years
    aux = t2m_anomaly_mean(:,:,i);
    
    for j = 1:5
        eval(['aux_cls',num2str(j),' = aux(CLASS_th==j);']);
        eval(['AUX_',num2str(j),' = [AUX_',num2str(j),',aux_cls',num2str(j),'];'])
    end
end 

bidedges = [-5:0.05:5];
f = figure;
f.Position = [100 100 800 600];
vals = cell(2,5);
for j=1:5
    subplot(5,4,(j-1).*4 +1)
    %     eval(['h = histogram(AUX_',num2str(j),'(:),''BinEdges'',bidedges,''Normalization'',''probability'',''DisplayStyle'',''stairs'',''EdgeColor'',colors_th(j,:)),''Linewidth'',5'])
    eval(['h = histogram(AUX_',num2str(j),'(:),''BinEdges'',bidedges,''Normalization'',''pdf'',''FaceColor'',colors_th(j,:),''EdgeColor'',colors_th(j,:))'])
    hold on 
    % Fit a normal distribution to the data
    eval(['data = AUX_',num2str(j),'(:);'])
    % Get the parameters of the fitted distribution
    mu = mean(data,'omitnan');
    sigma = std(data,'omitnan');
    set(gca,'Fontsize',18)
    grid on
    xlabel('Air T anomaly (\circC)')
    ylabel('pdf [-]')
    
    % Create a range of x values for plotting the fitted Gaussian
    x = bidedges;
    
    % Compute the PDF of the fitted Gaussian
    y = (1 / (sigma * sqrt(2 * pi))) * exp(-0.5 * ((x - mu) / sigma) .^ 2);
    
    MU(j)    = mu;
    SIGMA(j) = sigma;
    
    % Plot the fitted Gaussian on top of the histogram
    subplot(5,4,sort([[2,3,4,5].*4-1,[2,3,4,5].*4]))
    plot(x, y, 'color',colors_th(j,:), 'LineWidth', 4);
    hold on;
    
    H(j) = h;
    hold on
    set(gca,'Fontsize',18)
    grid on
    xlabel('Air Temperature anomaly (\circC)')
    ylabel('pdf [-]')
    title({'Air T anomaly distribution per thermokarst';'coverage classification'})
    
    vals(1,j) = num2cell(mu);
    vals(2,j) = num2cell(sigma);
end



% Create column names
columnNames = tickLabels_th;

% Create the table in the figure
uitable('Data', vals, 'ColumnName', columnNames, 'RowName', {'Mean', 'STD'}, ...
    'Position', [20 20 600 100],'FontSize', 18);

legend(tickLabels_th)
 
if saveif save_flag
    saveas(gcf,[path_fig,'Thermokarst_airTAnom_hist_subplot_table.png'])
    saveas(gcf,[path_fig,'Thermokarst_airTAnom_hist_subplot_table.fig'])
end 