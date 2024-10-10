% ###################################################################
% ####### plotPolarNDVI_EVI  ####
% This script reads the EVI and NDVI average values in June/july for 
% each year generated using the MODIS Terra and Aqua
% data and generate the polar plots. 

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
save_flag = 0;

% Read the NetCDF file 
filename = [path_in,'NDVI_EVI_yearly.nc'];
ncdisp(filename)

NDVI      = ncread(filename,'NDVI_Arctic');
% Filter negative and non valid values
NDVI(NDVI<=0) = NaN;

EVI       = ncread(filename,'EVI_Arctic');
% Filter negative and non valid values
EVI(EVI<=0) = NaN;

YEAR      = ncread(filename,'YEAR');
latitude  = ncread(filename,'latitude');
longitude = ncread(filename,'longitude');


if plotting == 1
    % Plotting NDVI per year (June-July) every 8-days
    f = figure;
    f.Position = [100 100 1800 2000];
    for i = 1:length(YEAR)
        subplot(4,4,i)
        m_proj('stereographic', 'lat',90,'long',30,'radius',30);
        aux = squeeze(NDVI(i,:,:));
        aux(aux<=0)=NaN; % Non valid         
        m_pcolor(longitude,latitude,aux');
        shading flat; % Optional: removes grid lines in pcolor plots
        m_coast('LineStyle','-','color','k');
        m_grid('box','fancy','tickdir','in','fontsize', 18);
        h = colorbar;
        % Adjust colorbar position
        cb_pos    = get(h, 'Position');
        cb_pos(1) = cb_pos(1) + 0.02; % Adjust the amount of shift as needed
        set(h, 'Position', cb_pos,'Fontsize',18)
        ylabel(h,'NDVI','Fontsize',18)
        t = title(['Year ', num2str(YEAR(i))]);
        set(gca,'Fontsize',18)
        set(t,'Position',get(t,'Position')+[0 .04 0]);  % move up slightly
        %         caxis([0.6,1])
    end
    if save_flag
        saveas(gcf,[path_fig,'NDVI_TerraAqua_JunJul_025deg_2003_2016.png'])
        saveas(gcf,[path_fig,'NDVI_TerraAqua_JunJul_025deg_2003_2016.fig'])
    end
    
    
    
    
    % Plotting EVI per year (June-July) every 8-days
    f = figure;
    f.Position = [100 100 1800 2000];
    for i = 1:length(YEAR)
        subplot(4,4,i)
        m_proj('stereographic', 'lat',90,'long',30,'radius',30);
        m_pcolor(longitude,latitude,squeeze(EVI(i,:,:))');
        shading flat; % Optional: removes grid lines in pcolor plots
        m_coast('LineStyle','-','color','k');
        m_grid('box','fancy','tickdir','in','fontsize', 18);
        h = colorbar;
        % Adjust colorbar position
        cb_pos    = get(h, 'Position');
        cb_pos(1) = cb_pos(1) + 0.02; % Adjust the amount of shift as needed
        set(h, 'Position', cb_pos,'Fontsize',18)
        ylabel(h,'EVI','Fontsize',18)
        t = title(['Year ', num2str(YEAR(i))]);
        set(gca,'Fontsize',18)
        set(t,'Position',get(t,'Position')+[0 .04 0]);  % move up slightly
        caxis([-0.2,0.6])
    end
    if save_flag
        saveas(gcf,[path_fig,'EVI_TerraAqua_JunJul_025deg_2003_2016.png'])
        saveas(gcf,[path_fig,'EVI_TerraAqua_JunJul_025deg_2003_2016.fig'])
    end
    
    
    %Transition between years
    cmap = [255/255, 0, 0; 151/255, 208/255, 119/255];
    f = figure;
    f.Position = [100 100 1800 2000];
    for i = 1:length(YEAR)-1
        subplot(4,4,i)
        m_proj('stereographic', 'lat',90,'long',30,'radius',30);
        aux = squeeze(EVI(i+1,:,:)-EVI(i,:,:));
        aux(aux<0) = -1;
        aux(aux>0) = 1;
        m_pcolor(longitude,latitude,aux');
        shading flat; % Optional: removes grid lines in pcolor plots
        m_coast('LineStyle','-','color','k');
        m_grid('box','fancy','tickdir','in','fontsize', 18);
        colormap(cmap)
        h = colorbar;
        % Adjust colorbar position
        cb_pos    = get(h, 'Position');
        cb_pos(1) = cb_pos(1) + 0.02; % Adjust the amount of shift as needed
        set(h, 'Position', cb_pos,'Fontsize',18)
        %ylabel(h,'\Delta EVI','Fontsize',18)        
        h.Ticks = [-0.5, 0.5]; % Positions of ticks
        h.TickLabels = {'EVI \downarrow ',' EVI \uparrow'}; % Custom labels
        h.FontSize = 16;        
        t = title(['Transition ', num2str(YEAR(i)),'-',num2str(YEAR(i+1))]);
        set(gca,'Fontsize',18)
        set(t,'Position',get(t,'Position')+[0 .04 0]);  % move up slightly
    end
    if save_flag
        saveas(gcf,[path_fig,'EVI_TerraAqua_JunJul_025deg_2003_2016_Transition.png'])
        saveas(gcf,[path_fig,'EVI_TerraAqua_JunJul_025deg_2003_2016_Transition.fig'])
    end

    
    %Transition between years
    cmap = [255/255, 0, 0; 151/255, 208/255, 119/255];
    f = figure;
    f.Position = [100 100 1800 2000];
    for i = 1:length(YEAR)-1
        subplot(4,4,i)        
        m_proj('stereographic', 'lat',90,'long',30,'radius',30);
        aux = squeeze(NDVI(i+1,:,:)-NDVI(i,:,:));
        aux(aux<0) = -1;
        aux(aux>0) = 1;
        m_pcolor(longitude,latitude,aux');
        shading flat; % Optional: removes grid lines in pcolor plots
        m_coast('LineStyle','-','color','k');
        m_grid('box','fancy','tickdir','in','fontsize', 18);
        colormap(cmap)
        h = colorbar;
        % Adjust colorbar position
        cb_pos    = get(h, 'Position');
        cb_pos(1) = cb_pos(1) + 0.02; % Adjust the amount of shift as needed
        set(h, 'Position', cb_pos,'Fontsize',18)
        %ylabel(h,'\Delta EVI','Fontsize',18)        
        h.Ticks = [-0.5, 0.5]; % Positions of ticks
        h.TickLabels = {'NDVI \downarrow ',' NDVI \uparrow'}; % Custom labels
        h.FontSize = 16;        
        t = title(['Transition ', num2str(YEAR(i)),'-',num2str(YEAR(i+1))]);
        set(gca,'Fontsize',18)
        set(t,'Position',get(t,'Position')+[0 .04 0]);  % move up slightly
    end
    if save_flag
        saveas(gcf,[path_fig,'NDVI_TerraAqua_JunJul_025deg_2003_2016_Transition.png'])
        saveas(gcf,[path_fig,'NDVI_TerraAqua_JunJul_025deg_2003_2016_Transition.fig'])
    end

    
    
    
    %Transition between years 2015-2016-2017 
    cmap = [255/255, 0, 0; 151/255, 208/255, 119/255];
    f = figure;
    f.Position = [100 100 1800 600];
    for i = 1:2
        subplot(1,2,i)
        m_proj('stereographic', 'lat',90,'long',30,'radius',30);
        aux = squeeze(NDVI(12+i+1,:,:)-NDVI(12+i,:,:));
        aux(aux<0) = -1;
        aux(aux>0) = 1;
        m_pcolor(longitude,latitude,aux');
        shading flat; % Optional: removes grid lines in pcolor plots
        m_coast('LineStyle','-','color','k');
        m_grid('box','fancy','tickdir','in','fontsize', 18);
        colormap(cmap)
        h = colorbar;
        % Adjust colorbar position
        cb_pos    = get(h, 'Position');
        cb_pos(1) = cb_pos(1) + 0.02; % Adjust the amount of shift as needed
        set(h, 'Position', cb_pos,'Fontsize',18)     
        h.Ticks = [-0.5, 0.5]; % Positions of ticks
        h.TickLabels = {'NDVI \downarrow ',' NDVI \uparrow'}; % Custom labels
        h.FontSize = 16;        
        t = title(['Transition ', num2str(YEAR(i+12)),'-',num2str(YEAR(i+12+1))]);
        set(gca,'Fontsize',18)
        set(t,'Position',get(t,'Position')+[0 .04 0]);  % move up slightly
    end
    
    
     %Transition between years 2015-2016-2017 
    cmap = [255/255, 0, 0; 151/255, 208/255, 119/255];
    f = figure;
    f.Position = [100 100 1200 500];
    for i = 1:2
        subplot(1,2,i)
        m_proj('stereographic', 'lat',90,'long',30,'radius',30);
        aux = squeeze(EVI(12+i+1,:,:)-EVI(12+i,:,:));
        aux(aux<0) = -1;
        aux(aux>0) = 1;
        m_pcolor(longitude,latitude,aux');
        shading flat; % Optional: removes grid lines in pcolor plots
        m_coast('LineStyle','-','color','k');
        m_grid('box','fancy','tickdir','in','fontsize', 18);
        colormap(cmap)
        h = colorbar;
        % Adjust colorbar position
        cb_pos    = get(h, 'Position');
        cb_pos(1) = cb_pos(1) + 0.07; % Adjust the amount of shift as needed
        set(h, 'Position', cb_pos,'Fontsize',18)
        %ylabel(h,'\Delta EVI','Fontsize',18)        
        h.Ticks = [-0.5, 0.5]; % Positions of ticks
        h.TickLabels = {'EVI \downarrow ',' NDVI \uparrow'}; % Custom labels
        h.FontSize = 16;        
        t = title(['Transition ', num2str(YEAR(i+12)),'-',num2str(YEAR(i+12+1))]);
        set(gca,'Fontsize',18)
        set(t,'Position',get(t,'Position')+[0 .04 0]);  % move up slightly
    end
    
    
    %Transition between years 2015-2016-2017 
    f = figure;
    f.Position = [100 100 1300 600];
    cmap = redblue(length([-0.1:0.01:0.1])); % EVI delta range
    cmap = [0.85 0.85 0.85; cmap];
    for i = 1:2
        subplot(1,2,i)
        m_proj('stereographic', 'lat',90,'long',30,'radius',30);
        aux = squeeze(EVI(12+i+1,:,:)-EVI(12+i,:,:));
        aux(isnan(aux))=-300;
        m_pcolor(longitude,latitude,aux');
        shading flat; % Optional: removes grid lines in pcolor plots
        m_coast('LineStyle','-','color','k');
        m_grid('box','fancy','tickdir','in','fontsize', 18);
        colormap(cmap)
        
        h = colorbar;
        % Adjust colorbar position
        cb_pos    = get(h, 'Position');
        cb_pos(1) = cb_pos(1) + 0.06; % Adjust the amount of shift as needed
        set(h, 'Position', cb_pos,'Fontsize',18)
%       ylabel(h,'\Delta NDVI','Fontsize',18)         
        
        % Set the title for the colorbar
        titleText = '\Delta EVI';       
        
        % Create a text object for the title above the colorbar
        titleHandle = title(h, titleText);
        
        % Adjust the position of the title to place it above the colorbar
        set(titleHandle, 'Units', 'normalized');
        titlePosition = get(titleHandle, 'Position');
        titlePosition(2) =  1.05; % Adjust this value to move the title up
        set(titleHandle, 'Position', titlePosition);
        
        % Adjust the font size and other properties as needed
        set(titleHandle, 'FontSize', 18, 'FontWeight', 'bold');
        
        set(get(h, 'Label'), 'String', 'NaN');
        lh = get(h, 'Label');
        labelPos = get(lh, 'Position');
        set(lh, 'Position', [1,-0.11], 'Rotation', 0, 'HorizontalAlignment', 'center');

        t = title(['Transition ', num2str(YEAR(i+12)),'-',num2str(YEAR(i+12+1))]);
        set(gca,'Fontsize',18)
        set(t,'Position',get(t,'Position')+[0 .04 0]);  % move up slightly
        caxis([-0.11,0.1])
        nanColorIndex = 1; % The index of the gray color in the new colormap
        set(gca, 'Color', cmap(nanColorIndex, :));
    end
    if save_flag
        saveas(gcf,[path_fig,'EVI_TerraAqua_JunJul_025deg_2016_Transitions_abs.png'])
        saveas(gcf,[path_fig,'EVI_TerraAqua_JunJul_025deg_2016_Transitions_abs.fig'])
    end
    
    
    %Transition between years 2015-2016-2017 
    f = figure;
    f.Position = [100 100 1400 600];
    cmap = redblue(length([-0.1:0.01:0.1])); % NDVI delta range
    cmap = [0.85 0.85 0.85; cmap];
%     cmap(find(sum(cmap,2)==3),:) = [0.5,0.5,0.5];% Replace 0 white per black
    for i = 1:2
        subplot(1,2,i)
        m_proj('stereographic', 'lat',90,'long',30,'radius',30);
        aux = squeeze(NDVI(12+i+1,:,:)-NDVI(12+i,:,:));
        aux(isnan(aux))=-300;
        m_pcolor(longitude,latitude,aux');
        shading flat; % Optional: removes grid lines in pcolor plots
        m_coast('LineStyle','-','color','k');
        m_grid('box','fancy','tickdir','in','fontsize', 18);
        colormap(cmap)
        
        h = colorbar;
        % Adjust colorbar position
        cb_pos    = get(h, 'Position');
        cb_pos(1) = cb_pos(1) + 0.06; % Adjust the amount of shift as needed
        set(h, 'Position', cb_pos,'Fontsize',18)
%       ylabel(h,'\Delta NDVI','Fontsize',18)         
        
        % Set the title for the colorbar
        titleText = '\Delta NDVI';       
        
        % Create a text object for the title above the colorbar
        titleHandle = title(h, titleText);
        
        % Adjust the position of the title to place it above the colorbar
        set(titleHandle, 'Units', 'normalized');
        titlePosition = get(titleHandle, 'Position');
        titlePosition(2) =  1.05; % Adjust this value to move the title up
        set(titleHandle, 'Position', titlePosition);
        
        % Adjust the font size and other properties as needed
        set(titleHandle, 'FontSize', 18, 'FontWeight', 'bold');
        
        set(get(h, 'Label'), 'String', 'NaN');
        lh = get(h, 'Label');
        labelPos = get(lh, 'Position');
        set(lh, 'Position', [1,-0.11], 'Rotation', 0, 'HorizontalAlignment', 'center');

        t = title(['Transition ', num2str(YEAR(i+12)),'-',num2str(YEAR(i+12+1))]);
        set(gca,'Fontsize',18)
        set(t,'Position',get(t,'Position')+[0 .04 0]);  % move up slightly
        caxis([-0.11,0.1])
        nanColorIndex = 1; % The index of the gray color in the new colormap
        set(gca, 'Color', cmap(nanColorIndex, :));
    end
    if save_flag
        saveas(gcf,[path_fig,'NDVI_TerraAqua_JunJul_025deg_2016_Transitions_abs.png'])
        saveas(gcf,[path_fig,'NDVI_TerraAqua_JunJul_025deg_2016_Transitions_abs.fig'])
    end
    
    %Transition between years 2015-2016-2017 
    f = figure;
    f.Position = [100 100 2300 600];
    for i = 1:3
        subplot(1,3,i)
        m_proj('stereographic', 'lat',90,'long',30,'radius',30);
        aux = squeeze(NDVI(12+i,:,:));
        m_pcolor(longitude,latitude,aux');
        shading flat; % Optional: removes grid lines in pcolor plots
        m_coast('LineStyle','-','color','k');
        m_grid('box','fancy','tickdir','in','fontsize', 18);
        h = colorbar;
        % Adjust colorbar position
        cb_pos    = get(h, 'Position');
        cb_pos(1) = cb_pos(1) + 0.04; % Adjust the amount of shift as needed
        set(h, 'Position', cb_pos,'Fontsize',18)
        ylabel(h,'\ NDVI','Fontsize',18)         
        t = title(['Year ', num2str(YEAR(i+12))]);
        set(gca,'Fontsize',18)
        set(t,'Position',get(t,'Position')+[0 .04 0]);  % move up slightly
        caxis([0.4,0.8])
    end
    
    
end






