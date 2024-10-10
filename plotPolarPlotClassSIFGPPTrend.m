% ###################################################################
% ####### PLOT POLAT PLOTS CLASSIFICATION SIF AND GPP TRENDS ####
% This script generates polar plots for the SIF & GPP classification trends:
% There are four classes combining the negative and positive GPP and SIF
% trends.
% Class  2 both +: GPP+, SIF +
% Class  0 gpp - : GPP-, SIF +
% Class  1 gpp+  :GPP+, SIF -
% Class -1 both -: GPP+, SIF +
% 
% %
% Author: Neus Sabater
% e-mail: neus.sabater@fmi.fi 
% version 0: (19/02/2024) 
% ------------------------------------------------------
% ------------------------------------------------------
% setting for saving and ploting 

plotting    = 1;
save_flag   = 0;
years_      = 'final';

addpath(genpath([pwd,filesep,'functions_environment']))
% Automatically detect environment
env = detect_environment();
% Load configuration
config = load_conf(env);

path_fig    = config.output_dir.general;
path_in     = config.input_dir.computed;




% ------------------------------------------------------
% Read the SIF GPP separate trend for the Arctic region and crop the region
% of Finland.
% ------------------------------------------------------
% ------------------------------------------------------
% read the slope trend values 
switch years_
    case 'initial'
      load([path_in,'SIF_GPP_Z_score_TrendArticSIFGPPFluxSat_JJ_2002_2016.mat'])   
    case 'final'
     load([path_in,'SIF_GPP_Z_score_TrendArticSIFGPPFluxSat_JJ_2003_2018.mat'])   
end 

% This file contains: 
% SIF_Trend
% GPP_Trend
% SIF_TrendStore
% GPP_TrendStore
% longitude
% latitude_range

[LAT,LON]                = meshgrid(latitude_range,longitude); % coordinates of the SIF and GPP data 
SIFTrend_evFIN           = squeeze(SIF_Trend(:,:,2))';
GPPTrend_evFIN           = squeeze(GPP_Trend(:,:,2))';
SIFTrend_Zscore_evFIN    = squeeze(SIF_Zscore_Trend(:,:,2))';
GPPTrend_Zscore_evFIN    = squeeze(GPP_Zscore_Trend(:,:,2))';


classes  = NaN(size(GPPTrend_evFIN));
classes(GPPTrend_evFIN>=0 & SIFTrend_evFIN>=0)=2;
classes(GPPTrend_evFIN<=0 & SIFTrend_evFIN<=0)=-1;
classes(GPPTrend_evFIN>=0 & SIFTrend_evFIN<=0)=1;
classes(GPPTrend_evFIN<=0 & SIFTrend_evFIN>=0)=0; 


%% mapping in polar plots
% ---------------------------
addpath(genpath(config.matlab_path.mappingToolbox))
    
if plotting == 1
    f = figure;
    f.Position = [100 100 1500 1000];
    m_proj('stereographic', 'lat',90,'long',30,'radius',30); 
    m_pcolor(longitude,latitude_range,classes); % Ensure lon and lat match the grid of land_cover data
    shading flat; % Optional: removes grid lines in pcolor plots
    m_coast('LineStyle','-','color','k');
    m_grid('box','fancy','tickdir','in','fontsize', 18);
    h = colorbar;
    h.Ticks = [-1, 0, 1,2]; % Positions of ticks
    h.TickLabels = {'GPP\downarrow - SIF\downarrow ', 'GPP\downarrow - SIF\uparrow', 'GPP\uparrow - SIF\downarrow','GPP\uparrow - SIF\uparrow'}; % Custom labels
    switch years_
        case 'initial'
            ylabel(h,'GPP - SIF trend 2002-2016')
            set(h,'Fontsize',18)
            t = title('GPP - SIF trend 2002-2016');
        case 'final'
            ylabel(h,'GPP - SIF trend 2003-2018')
            set(h,'Fontsize',18)
            t = title('GPP - SIF trend 2003-2018');
    end
    set(t,'Position',get(t,'Position')+[0 .03 0],'Fontsize',18);  % move up slightly
    caxis([-1.5,2.5])
    cmap = [255/255, 0, 0; 0,0,0; 1,1,0.5; 151/255, 208/255, 119/255];
    colormap(cmap)
    
    % Adjust colorbar position
    cb_pos = get(h, 'Position');
    cb_pos(1) = cb_pos(1) + 0.08; % Adjust the amount of shift as needed
    set(h, 'Position', cb_pos);
    
    if save_flag
        switch years_
            case 'initial'
                saveas(gcf,[path_fig,'GPPSIFratiotrend_Arctic_005deg_2002_2016.png'])
                saveas(gcf,[path_fig,'GPPSIFratiotrend_Arctic_005deg_2002_2016.fig'])
                
            case 'final'
                saveas(gcf,[path_fig,'GPPSIFratiotrend_Arctic_005deg_2003_2018.png'])
                saveas(gcf,[path_fig,'GPPSIFratiotrend_Arctic_005deg_2003_2018.fig'])
        end
    end
end

    



%% Compute the anomalies per year
% ---------------------------
SIF_mean           = mean(SIF_TrendStore,3,'omitnan');
SIF_median         = median(SIF_TrendStore,3,'omitnan');
SIF_anomaly_mean   = SIF_TrendStore-SIF_mean;
SIF_anomaly_median = SIF_TrendStore-SIF_median;


GPP_mean           = mean(GPP_TrendStore,3,'omitnan');
GPP_median         = median(GPP_TrendStore,3,'omitnan');
GPP_anomaly_mean   = GPP_TrendStore-GPP_mean;
GPP_anomaly_median = GPP_TrendStore-GPP_median;

if save_flag
    switch years_
        case 'initial'
            years_covered_GPP = [2002:2016];
            years_covered_SIF = [2002:2016];
            save('GPPAnomaly_JunJul_005deg_2002_2016.mat','GPP_anomaly_mean','GPP_anomaly_median',...
                'GPP_median','GPP_mean','longitude','latitude_range','years_covered_GPP')
            save('SIFAnomaly_JunJul_005deg_2002_2016.mat','SIF_anomaly_mean','SIF_anomaly_median',...
                'SIF_median','SIF_mean','longitude','latitude_range','years_covered_SIF')
            
        case 'final'
            years_covered_GPP = [2003:2018];
            years_covered_SIF = [2003:2018];
            save('GPPAnomaly_JunJul_005deg_2003_2018.mat','GPP_anomaly_mean','GPP_anomaly_median',...
                'GPP_median','GPP_mean','longitude','latitude_range','years_covered_GPP')
            save('SIFAnomaly_JunJul_005deg_2003_2018.mat','SIF_anomaly_mean','SIF_anomaly_median',...
                'SIF_median','SIF_mean','longitude','latitude_range','years_covered_GPP')
    end    
end


%% mapping in polar plots of anomalies 
% ---------------------------
    
if plotting == 1
    f = figure;
    f.Position = [100 100 2000 2000];
    for i=1:size(GPP_TrendStore,3)
        subplot(4,4,i)
        m_proj('stereographic', 'lat',90,'long',30,'radius',30);
        m_pcolor(longitude,latitude_range,squeeze(SIF_anomaly_mean(:,:,i))'); % Ensure lon and lat match the grid of land_cover data
        shading flat; % Optional: removes grid lines in pcolor plots
        m_coast('LineStyle','-','color','k');
        m_grid('box','fancy','tickdir','in','fontsize', 18);
        h         = colorbar;
        caxis([-0.1,0.1])
        % Adjust colorbar position
        cb_pos    = get(h, 'Position');
        cb_pos(1) = cb_pos(1) + 0.02; % Adjust the amount of shift as needed
        set(h, 'Position', cb_pos,'Fontsize',18)
        ylabel(h,'SIF Anomaly','Fontsize',18)
        cmap        = redblue(length(-0.1:0.02:0.1));
        [row,~]     = find(sum(cmap,2)==3);
        cmap(row,:) = [0.6,0.6,0.6];
        colormap(cmap);
        switch years_
            case 'initial'
                title([num2str(i+2001)],'Fontsize',18) % Started in 2002
            case 'final'
                title([num2str(i+2002)],'Fontsize',18) % Started in 2003
        end
    end
    
    if save_flag
        switch years_
            case 'initial'
                saveas(gcf,[path_fig,'SIFanomaly_Arctic_005deg_2002_2016.png'])
                saveas(gcf,[path_fig,'SIFanomaly_Arctic_005deg_2002_2016.fig'])
                
            case 'final'
                saveas(gcf,[path_fig,'SIFanomaly_Arctic_005deg_2003_2018.png'])
                saveas(gcf,[path_fig,'SIFanomaly_Arctic_005deg_2003_2018.fig'])
        end
    end
end

%% GPP anomalies
if plotting == 1
    f = figure;
    f.Position = [100 100 2000 2000];
    for i=1:size(GPP_TrendStore,3)
        subplot(4,4,i)
        m_proj('stereographic', 'lat',90,'long',30,'radius',30);
        m_pcolor(longitude,latitude_range,squeeze(GPP_anomaly_mean(:,:,i))'); % Ensure lon and lat match the grid of land_cover data
        shading flat; % Optional: removes grid lines in pcolor plots
        m_coast('LineStyle','-','color','k');
        m_grid('box','fancy','tickdir','in','fontsize', 18);
        h         = colorbar;
        caxis([-4,4])
        % Adjust colorbar position
        cb_pos    = get(h, 'Position');
        cb_pos(1) = cb_pos(1) + 0.02; % Adjust the amount of shift as needed
        set(h, 'Position', cb_pos,'Fontsize',18)
        ylabel(h,'GPP Anomaly','Fontsize',18)
        cmap        = redblue(length(-0.1:0.02:0.1));
        [row,~]     = find(sum(cmap,2)==3);
        cmap(row,:) = [0.6,0.6,0.6];
        colormap(cmap);
        switch years_
            case 'initial'
                title([num2str(i+2001)],'Fontsize',18) % Started in 2002
            case 'final'
                title([num2str(i+2002)],'Fontsize',18) % Started in 2003
        end
    end
    
    if save_flag
        switch years_
            case 'initial'
                saveas(gcf,[path_fig,'GPPanomaly_Arctic_005deg_2002_2016.png'])
                saveas(gcf,[path_fig,'GPPanomaly_Arctic_005deg_2002_2016.fig'])
                
            case 'final'
                saveas(gcf,[path_fig,'GPPanomaly_Arctic_005deg_2003_2018.png'])
                saveas(gcf,[path_fig,'GPPanomaly_Arctic_005deg_2003_2018.fig'])
        end
    end
end





%% Compute the GPP SIF classification per year
% ---------------------------
if plotting == 1
    f = figure;
    f.Position = [100 100 2000 1200];
    
    for i= 1:size(GPP_TrendStore,3)-1 % num of years
        
        SIFTrend_evFIN           = (abs(squeeze(SIF_TrendStore(:,:,i+1)))-abs(squeeze(SIF_TrendStore(:,:,i))))';
        GPPTrend_evFIN           = (abs(squeeze(GPP_TrendStore(:,:,i+1)))-abs(squeeze(GPP_TrendStore(:,:,i))))';
        
        
        classes  = NaN(size(GPPTrend_evFIN));
        classes(GPPTrend_evFIN>=0 & SIFTrend_evFIN>=0)=2;
        classes(GPPTrend_evFIN<=0 & SIFTrend_evFIN<=0)=-1;
        classes(GPPTrend_evFIN>=0 & SIFTrend_evFIN<=0)=1;
        classes(GPPTrend_evFIN<=0 & SIFTrend_evFIN>=0)=0;
        
        
        
        subplot(4,4,i)
        m_proj('stereographic', 'lat',90,'long',30,'radius',30);
        m_pcolor(longitude,latitude_range,classes); % Ensure lon and lat match the grid of land_cover data
        shading flat; % Optional: removes grid lines in pcolor plots
        m_coast('LineStyle','-','color','k');
        m_grid('box','fancy','tickdir','in','fontsize', 18);
        h = colorbar;
        h.Ticks = [-1, 0, 1,2]; % Positions of ticks
        h.TickLabels = {'GPP\downarrow - SIF\downarrow ', 'GPP\downarrow - SIF\uparrow', 'GPP\uparrow - SIF\downarrow','GPP\uparrow - SIF\uparrow'}; % Custom labels
        h.FontSize = 16;
        switch years_
            case 'initial'
                t =title(['\Delta: ',num2str(i+2001+1),'-',num2str(i+2001)],'Fontsize',18); % Started in 2002
            case 'final'
                t =title(['\Delta: ',num2str(i+2002+1),'-',num2str(i+2002)],'Fontsize',18); % Started in 2003
        end
        set(t,'Position',get(t,'Position')+[0 .1 0],'Fontsize',18);  % move up slightly
        caxis([-1.5,2.5])
        cmap = [255/255, 0, 0; 0,0,0; 1,1,0.5; 151/255, 208/255, 119/255];
        colormap(cmap)
        
        %     Adjust colorbar position
        cb_pos = get(h, 'Position');
        cb_pos(1) = cb_pos(1) + 0.04; % Adjust the amount of shift as needed
        set(h, 'Position', cb_pos);
end
    
    
  if save_flag
        switch years_
            case 'initial'
                saveas(gcf,[path_fig,'SIFGPPyearly_classification_2002_2016.png'])
                saveas(gcf,[path_fig,'SIFGPPyearly_classification_2002_2016.fig'])
                
            case 'final'
                saveas(gcf,[path_fig,'SIFGPPyearly_classification_2003_2018.png'])
                saveas(gcf,[path_fig,'SIFGPPyearly_classification_2003_2018.fig'])
        end
    end


end 



%% Compute the GPP SIF classification per year
% and the  statistics of the underlaying cover
% ---------------------------
BOR_int = ncread([path_in,'BAWLD_dataset.nc'],'BOR_int');
GLA_int = ncread([path_in,'BAWLD_dataset.nc'],'GLA_int');
TUN_int = ncread([path_in,'BAWLD_dataset.nc'],'TUN_int');
ROC_int = ncread([path_in,'BAWLD_dataset.nc'],'ROC_int');
WET_int = ncread([path_in,'BAWLD_dataset.nc'],'WET_int');
LAK_int = ncread([path_in,'BAWLD_dataset.nc'],'LAK_int');
RIV_int = ncread([path_in,'BAWLD_dataset.nc'],'RIV_int');
CLASS_int = ncread([path_in,'BAWLD_dataset.nc'],'CLASS_int');




BOR_stat_c1 = []; BOR_stat_cminus1 = []; BOR_stat_c2 = []; BOR_stat_c0 = [];
WET_stat_c1 = []; WET_stat_cminus1 = []; WET_stat_c2 = []; WET_stat_c0 = [];
TUN_stat_c1 = []; TUN_stat_cminus1 = []; TUN_stat_c2 = []; TUN_stat_c0 = [];
CLASS_stat_c1 = []; CLASS_stat_cminus1 = []; CLASS_stat_c2 = []; CLASS_stat_c0 = [];


for i= 1:size(GPP_TrendStore,3)-1 % num of years
    
    SIFTrend_evFIN           = (abs(squeeze(SIF_TrendStore(:,:,i+1)))-abs(squeeze(SIF_TrendStore(:,:,i))));
    GPPTrend_evFIN           = (abs(squeeze(GPP_TrendStore(:,:,i+1)))-abs(squeeze(GPP_TrendStore(:,:,i))));
    
    
    classes  = NaN(size(BOR_int));
    classes(GPPTrend_evFIN>=0 & SIFTrend_evFIN>=0)=2;
    classes(GPPTrend_evFIN<=0 & SIFTrend_evFIN<=0)=-1;
    classes(GPPTrend_evFIN>=0 & SIFTrend_evFIN<=0)=1;
    classes(GPPTrend_evFIN<=0 & SIFTrend_evFIN>=0)=0;
    
     
    % BOREAL 
    BOR_stat_c1 = [BOR_stat_c1,BOR_int(classes==1)'];
    BOR_stat_c2 = [BOR_stat_c2,BOR_int(classes==2)'];
    BOR_stat_cminus1 = [BOR_stat_cminus1,BOR_int(classes==-1)'];
    BOR_stat_c0 = [BOR_stat_c0,BOR_int(classes==0)'];
    
    % WETLAND
    WET_stat_c1 = [WET_stat_c1,WET_int(classes==1)'];
    WET_stat_c2 = [WET_stat_c2,WET_int(classes==2)'];
    WET_stat_cminus1 = [WET_stat_cminus1,WET_int(classes==-1)'];
    WET_stat_c0 = [WET_stat_c0,WET_int(classes==0)'];
    
    % TUNDRA
    TUN_stat_c1 = [TUN_stat_c1,TUN_int(classes==1)'];
    TUN_stat_c2 = [TUN_stat_c2,TUN_int(classes==2)'];
    TUN_stat_cminus1 = [TUN_stat_cminus1,TUN_int(classes==-1)'];
    TUN_stat_c0 = [TUN_stat_c0,TUN_int(classes==0)'];
    
    
    % Classes  k-means
    CLASS_stat_c1 = [CLASS_stat_c1,CLASS_int(classes==1)'];
    CLASS_stat_c2 = [CLASS_stat_c2,CLASS_int(classes==2)'];
    CLASS_stat_cminus1 = [CLASS_stat_cminus1,CLASS_int(classes==-1)'];
    CLASS_stat_c0 = [CLASS_stat_c0,CLASS_int(classes==0)'];
end



if plotting == 1
    
    f = figure;
    f.Position = [100 100 800 600];
    % Plot the stack bar plot for each class
    cmap = [255/255, 0, 0; 0,0,0; 1,1,0.5; 151/255, 208/255, 119/255];
    
    kmeans_c1 = histogram(CLASS_stat_c1(~isnan(CLASS_stat_c1)));
    kmeans_C1 = kmeans_c1.Values;
    kmeans_c2 = histogram(CLASS_stat_c2(~isnan(CLASS_stat_c2)));
    kmeans_C2 = kmeans_c2.Values;
    kmeans_cminus1 = histogram(CLASS_stat_cminus1(~isnan(CLASS_stat_cminus1)));
    kmeans_Cminus1 = kmeans_cminus1.Values;
    kmeans_c0 = histogram(CLASS_stat_c0(~isnan(CLASS_stat_c0)));
    kmeans_C0 = kmeans_c0.Values;
    
    total = sum([kmeans_C1;kmeans_C2;kmeans_Cminus1;kmeans_C0],1); % Normalize so the sum of each class is 1

    b = bar(1:15, (100.*[kmeans_C1;kmeans_C2;kmeans_Cminus1;kmeans_C0]./total)', 1, 'stack', 'FaceColor','flat','EdgeColor',[0.6,0.6,0.6]);
    ylim([0,100])
    b(1).CData = cmap(3,:);
    b(2).CData = cmap(4,:);
    b(3).CData = cmap(1,:);
    b(4).CData = cmap(2,:);
    legend('GPP\uparrow SIF \downarrow','GPP\uparrow SIF \uparrow ','GPP\downarrow SIF \downarrow','GPP\downarrow SIF \uparrow','Location','best')
    set(gca,'Fontsize',20)
    xticks([1:15])
    xticklabels({'Permafrost Peatlands','Sparse Boreal Peatlands','Rivers',...
        'Glaciers','Upland Tundra','Common Boreal Peatlands','Large Lakes',...
        'Lake-rich wetlands','Dominant Boreal Peatlands','Wetland-rich Tundra',...
        'Alpine and Tundra Barrens','Wetland and Lake-rich Tundra','Lake-rich Shield','Upland Boreal','Wetland and Lake-rich Yedoma Tundra'});
    xtickangle(45)
    ylabel('Percent of class presence (%)')
    grid off
    box on
    title('Arctic database distribution per GPP - SIF trend classification')
    
    
    if save_flag
        switch years_
            case 'initial'
                saveas(gcf,[path_fig,'BAWLD_SIFGPP_class_2002_2016.png'])
                saveas(gcf,[path_fig,'BAWLD_SIFGPP_class_2002_2016.fig'])
                
            case 'final'
                saveas(gcf,[path_fig,'BAWLD_SIFGPP_class_2003_2018.png'])
                saveas(gcf,[path_fig,'BAWLD_SIFGPP_class_2003_2018.fig'])
        end
    end
    
    
    figure 
    variableNames = {'Permafrost Peatlands','Sparse Boreal Peatlands','Rivers',...
        'Glaciers','Upland Tundra','Common Boreal Peatlands','Large Lakes',...
        'Lake-rich wetlands','Dominant Boreal Peatlands','Wetland-rich Tundra',...
        'Alpine and Tundra Barrens','Wetland and Lake-rich Tundra','Lake-rich Shield','Upland Boreal','Wetland and Lake-rich Yedoma Tundra'};
    classNames  = {'GPP\uparrow SIF \downarrow','GPP\uparrow SIF \uparrow ','GPP\downarrow SIF \downarrow','GPP\downarrow SIF \uparrow'};
    
    h = heatmap((100.*[kmeans_C1;kmeans_C2;kmeans_Cminus1;kmeans_C0]./total),'Title','Arctic database distribution per GPP - SIF trend classification',...
         'Colormap', bone);   
    % Set the custom x and y axis labels
    h.XDisplayLabels = variableNames;
    h.YDisplayLabels = classNames;
    h.CellLabelFormat = '%.0f';
    set(gca,'Fontsize',18)
    
    if save_flag
        switch years_
            case 'initial'
                saveas(gcf,[path_fig,'BAWLD_SIFGPP_class_2002_2016_heatmap.png'])
                saveas(gcf,[path_fig,'BAWLD_SIFGPP_class_2002_2016_heatmap.fig'])
                
            case 'final'
                saveas(gcf,[path_fig,'BAWLD_SIFGPP_class_2003_2018_heatmap.png'])
                saveas(gcf,[path_fig,'BAWLD_SIFGPP_class_2003_2018_heatmap.fig'])
        end
    end
    
end

%% Compute the GPP SIF anomaly per year
% for each of the underlaying cover
% ---------------------------
load([path_in,'airTAnomaly_JunJul_005deg_2003_2018.mat']) % 'years_cov','t2m_summer_anomaly_year','lon','lat'))

% if plotting == 1
    f = figure;
    f.Position = [100 100 2000 1200];
    
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
    
    % GPP anomalies scatter plots vs. Air T anomalies 
    
    for j=1:15
        auxT_an_store   = [];
        auxGPP_an_store = [];
        
        for i=1:size(GPP_TrendStore,3)
            auxGPP = squeeze(GPP_anomaly_mean(:,:,i));
            auxT   = squeeze(t2m_anomaly_mean(:,:,i));
            
            auxT   = auxT(CLASS_int==j);
            auxGPP = auxGPP(CLASS_int==j);
            
            auxGPP_an_store = [auxGPP_an_store(:);auxGPP(:)]; % all years
            auxT_an_store   = [auxT_an_store(:);auxT(:)]; % all years
        end
        subplot(4,4,j)
        d = scatter(auxT_an_store,auxGPP_an_store,8,colors(j,:),'filled');
        hold on
        d.MarkerFaceAlpha = 0.2;
        set(gca,'Fontsize',18)
        xlabel('Air T anomaly')
        ylabel('GPP anomaly')
        title(tickLabels{j})
        ylim([-5,5])
        xlim([-5,5])
        grid on
        box on
        x = auxT_an_store(~isnan(auxGPP_an_store) & ~isnan(auxT_an_store));
        y = auxGPP_an_store(~isnan(auxGPP_an_store) & ~isnan(auxT_an_store));
        X = [ones(length(x),1) x];
        b = X\y;
        GPP_B(j,:) = b;close all
        yCalc2 = X*b;
        plot(x,yCalc2,'k--')
        Rsq2 = 1 - sum((y - yCalc2).^2)/sum((y - mean(y)).^2);
        box on
        grid on
        grid minor
        set(gca,'TickLength',[0.05, 0.01])
        set(gca,'XMinorTick','on','YMinorTick','on')
        text(-4,3,['R^2 =',num2str(round(Rsq2,2))],'Fontsize',18)
        text(-4,4,['y = ',num2str(round(b(1),4)),' + ',num2str(round(b(2),4)),'x'],'Fontsize',18)
    end
    
    if save_flag
        switch years_
            case 'initial'
                saveas(gcf,[path_fig,'BAWLD_GPP_class_scatter_2002_2016.png'])
                saveas(gcf,[path_fig,'BAWLD_GPP_class_scatter_2002_2016.fig'])
                
            case 'final'
                saveas(gcf,[path_fig,'BAWLD_GPP_class_scatter_2003_2018.png'])
                saveas(gcf,[path_fig,'BAWLD_GPP_class_scatter_2003_2018.fig'])
        end
    end

    % SIF anomalies scatter plots vs. Air T anomalies 
    
    
    f = figure;
    f.Position = [100 100 2000 1200];
    for j=1:15
        auxT_an_store   = [];
        auxSIF_an_store = [];
        
        for i=1:size(GPP_TrendStore,3)
            auxSIF = squeeze(SIF_anomaly_mean(:,:,i));
            auxT   = squeeze(t2m_anomaly_mean(:,:,i));
            
            auxT   = auxT(CLASS_int==j);
            auxSIF = auxSIF(CLASS_int==j);
            
            auxSIF_an_store = [auxSIF_an_store(:);auxSIF(:)]; % all years
            auxT_an_store   = [auxT_an_store(:);auxT(:)]; % all years
        end
        subplot(4,4,j)
        d = scatter(auxT_an_store,auxSIF_an_store,8,colors(j,:),'filled');
        hold on
        d.MarkerFaceAlpha = 0.2;
        set(gca,'Fontsize',18)
        xlabel('Air T anomaly')
        ylabel('SIF anomaly')
        title(tickLabels{j})
        ylim([-0.3,0.3])
        xlim([-5,5])
        grid on
        box on
        x = auxT_an_store(~isnan(auxSIF_an_store) & ~isnan(auxT_an_store));
        y = auxSIF_an_store(~isnan(auxSIF_an_store) & ~isnan(auxT_an_store));
        X = [ones(length(x),1) x];
        b = X\y;
        SIF_B(j,:) = b;
        yCalc2 = X*b;
        plot(x,yCalc2,'k--')
        Rsq2 = 1 - sum((y - yCalc2).^2)/sum((y - mean(y)).^2);
        box on
        grid on
        grid minor
        set(gca,'TickLength',[0.05, 0.01])
        set(gca,'XMinorTick','on','YMinorTick','on')
        text(-4,0.1,['R^2 =',num2str(round(Rsq2,2))],'Fontsize',18)
        text(-4,0.2,['y = ',num2str(round(b(1),4)),' + ',num2str(round(b(2),4)),'x'],'Fontsize',18)
    end
    
    if save_flag
        switch years_
            case 'initial'
                saveas(gcf,[path_fig,'BAWLD_SIF_class_scatter_2002_2016.png'])
                saveas(gcf,[path_fig,'BAWLD_SIF_class_scatter_2002_2016.fig'])
                
            case 'final'
                saveas(gcf,[path_fig,'BAWLD_SIF_class_scatter_2003_2018.png'])
                saveas(gcf,[path_fig,'BAWLD_SIF_class_scatter_2003_2018.fig'])
        end
    end
    
% end

%% Figure of the anomalies rates of SIF and GPP per class
% from the linear regression

% if plotting == 1
    x = [-5:0.1:5]';
    X = [ones(length(x),1) x];
    for j=1:15
        f = figure(20);
        plot(x,X*SIF_B(j,:)','color',colors(j,:),'linewidth',3);
        hold on
        g = figure(40);
        plot(x,X*GPP_B(j,:)','color',colors(j,:),'linewidth',3);
        hold on
    end
    f = figure(20);
    f.Position = [100 100 1500 600];
    set(gca,'Fontsize',18)
    xlabel('Anomaly air T year^{-1}')
    ylabel('Anomaly SIF year^{-1}')
    box on
    grid on
    grid minor
    set(gca,'TickLength',[0.05, 0.01])
    set(gca,'XMinorTick','on','YMinorTick','on')
    legend(tickLabels,'Location','eastoutside')
    title('SIF anomalies vs. air Temperature anomaly')
    
    if save_flag
        switch years_
            case 'initial'
                saveas(gcf,[path_fig,'BAWLD_SIF_linearreg_2002_2016.png'])
                saveas(gcf,[path_fig,'BAWLD_SIF_linearreg_2002_2016.fig'])
                
            case 'final'
                saveas(gcf,[path_fig,'BAWLD_SIF_linearreg_2003_2018.png'])
                saveas(gcf,[path_fig,'BAWLD_SIF_linearreg_2003_2018.fig'])
        end
    end
    
    g = figure(40);
    g.Position = [100 100 1500 600];
    set(gca,'Fontsize',18)
    xlabel('Anomaly air T year^{-1}')
    ylabel('Anomaly GPP year^{-1}')
    box on
    grid on
    grid minor
    set(gca,'TickLength',[0.05, 0.01])
    set(gca,'XMinorTick','on','YMinorTick','on')
    legend(tickLabels,'Location','eastoutside')
    title('GPP anomalies vs. air Temperature anomaly')
    if save_flag
        switch years_
            case 'initial'
                saveas(gcf,[path_fig,'BAWLD_GPP_linearreg_2002_2016.png'])
                saveas(gcf,[path_fig,'BAWLD_GPP_linearreg_2002_2016.fig'])
                
            case 'final'
                saveas(gcf,[path_fig,'BAWLD_GPP_linearreg_2003_2018.png'])
                saveas(gcf,[path_fig,'BAWLD_GPP_linearreg_2003_2018.fig'])
        end
    end
    
% end



if plotting == 1
    
    figure
    for j=1:15
        yyaxis left
        h1 = plot(j,SIF_B(j,2),'o','Markersize',26,'markerfacecolor',colors(j,:),'markeredgecolor','k');
        hold on
        yyaxis right
        h2 = plot(j,GPP_B(j,2),'d','Markersize',26,'markerfacecolor',colors(j,:),'markeredgecolor','k');
        hold on
    end
    yyaxis left
    plot([1:15],SIF_B(:,2),'k-')
    yyaxis right
    plot([1:15],GPP_B(:,2),'k-')
    yyaxis right
    ylabel('GPP anomaly vs. air T anomaly')
    yyaxis left
    ylabel('SIF anomaly vs. air T anomaly')
    xticks([1:15])
    xticklabels({'Permafrost Peatlands','Sparse Boreal Peatlands','Rivers',...
        'Glaciers','Upland Tundra','Common Boreal Peatlands','Large Lakes',...
        'Lake-rich wetlands','Dominant Boreal Peatlands','Wetland-rich Tundra',...
        'Alpine and Tundra Barrens','Wetland and Lake-rich Tundra','Lake-rich Shield','Upland Boreal','Wetland and Lake-rich Yedoma Tundra'});
    xtickangle(45)
    box on
    grid on
    grid minor
    set(gca,'Fontsize',18)
    set(gca,'TickLength',[0.05, 0.005])
    set(gca,'XMinorTick','on','YMinorTick','on')
    xlim([0,16])
    title('GPP and SIF anomalies response to air Temperature')
    lgnd = legend([h1,h2],'SIF','GPP','Location','west','orientation','horizontal');
    set(lgnd,'color',[0.8,0.8,0.8]);
    if save_flag
        switch years_
            case 'initial'
                saveas(gcf,[path_fig,'BAWLD_GPPSIF_anom_changevsairT_2002_2016.png'])
                saveas(gcf,[path_fig,'BAWLD_GPPSIF_anom_changevsairT_2002_2016.fig'])
                
            case 'final'
                saveas(gcf,[path_fig,'BAWLD_GPPSIF_anom_changevsairT_2003_2018.png'])
                saveas(gcf,[path_fig,'BAWLD_GPPSIF_anom_changevsairT_2003_2018.fig'])
        end
    end
    
end


% Order the data in descend order
% if plotting == 1
    [~,posSIF] = sort(SIF_B(:,2),'descend');
    [~,posGPP]  = sort(GPP_B(:,2),'descend');
        
    
    tickLabels = {'Permafrost Peatlands','Sparse Boreal Peatlands','Rivers',...
        'Glaciers','Upland Tundra','Common Boreal Peatlands','Large Lakes',...
        'Lake-rich wetlands','Dominant Boreal Peatlands','Wetland-rich Tundra',...
        'Alpine and Tundra Barrens','Wetland and Lake-rich Tundra','Lake-rich Shield','Upland Boreal','Wetland and Lake-rich Yedoma Tundra'};
    
          for j=1:15
              tickLabels_ord{j} = tickLabels{posGPP(j)};
          end

          
          % Ordered using GPP
          f = figure;
          f.Position = [100 100 1500 600];
          for j=1:15
              yyaxis right
              %     h2 = plot(j,GPP_B(posGPP(j),2),'d','Markersize',30,'markerfacecolor',colors(posGPP(j),:),'markeredgecolor','k','Facealpha',0.5);
              h2 = bar(j,GPP_B(posGPP(j),2),'FaceAlpha',0.5,'Facecolor',colors(posGPP(j),:));
              
              hold on
              yyaxis left
              h1 = plot(j,SIF_B(posGPP(j),2),'o','Markersize',20,'markerfacecolor',colors(posGPP(j),:),'markeredgecolor','k');
              hold on
          end
          yyaxis right
          ylabel('GPP anomaly vs. air T anomaly')
          yyaxis left
          ylabel('SIF anomaly vs. air T anomaly')
          box on
          grid on
          ax = gca;
          ax.YAxis(1).Color = 'k';
          ax.YAxis(2).Color = 'k';
          xticks([1:15])
          xticklabels(tickLabels_ord)
          xtickangle(45)
          grid minor
          set(gca,'Fontsize',18)
          set(gca,'TickLength',[0.05, 0.005])
          set(gca,'XMinorTick','on','YMinorTick','on')
          xlim([0,16])
          title('GPP and SIF anomalies response to air Temperature')
          lgnd = legend([h1,h2],'SIF','GPP','Location','east','orientation','horizontal');
          set(lgnd,'color',[0.8,0.8,0.8]);
          
          if save_flag
              switch years_
                  case 'initial'
                      saveas(gcf,[path_fig,'BAWLD_GPPSIF_anom_changevsairT_order_2002_2016.png'])
                      saveas(gcf,[path_fig,'BAWLD_GPPSIF_anom_changevsairT_order_2002_2016.fig'])
                      
                  case 'final'
                      saveas(gcf,[path_fig,'BAWLD_GPPSIF_anom_changevsairT_order_2003_2018.png'])
                      saveas(gcf,[path_fig,'BAWLD_GPPSIF_anom_changevsairT_order_2003_2018.fig'])
              end
          end
% end

