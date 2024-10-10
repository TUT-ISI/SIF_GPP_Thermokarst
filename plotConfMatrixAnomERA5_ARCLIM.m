% #######################################################################################
% ####### PLOT CONFUSION MATRIX ANOMALIES OF SIF, GPP AND ERA5 AND ARCLIM VARIABLES  ####
% #######################################################################################
% This script generates confusion matrix between the anomalies of the
% different ERA5 and ARCLIM variables and the anomalies of GPP and SIF with the aim to 
% detect is the response is in sync.
% 
% %
% Author: Neus Sabater
% e-mail: neus.sabater@fmi.fi 
% version 0: (12/06/2024) 
% ------------------------------------------------------
% ------------------------------------------------------
% setting for saving and ploting 

addpath(genpath([pwd,filesep,'functions_environment']))
% Automatically detect environment
env = detect_environment();
% Load configuration
config = load_conf(env);

plotting    = 1;
save_flag   = 0;
years_      = 'final';
path_fig    = [config.output_dir.general,'ERA5_ARCLIM']; %'/Users/neus/Dropbox/IPL/02_PROYECTOS_PROPIOS/SIF_GPP/figures_all/figs/';
if ~exist(path_fig)
    mkdir(path_fig)
end

% ------------------------------------------------------
% Read the SIF GPP separate trend for the Arctic region 
% and crop the region of Finland.
% ------------------------------------------------------
% ------------------------------------------------------
path_in = config.input_dir.computed; % '/Users/neus/Dropbox/IPL/02_PROYECTOS_PROPIOS/SIF_GPP/codes/data/Computed_data/';

switch years_
    case 'initial'        
       disp('Not computed anomalies for this interval.')
       
    case 'final'
        filename = [path_in,'ERA5_ARCLIM_mean_anomalies_2003_2018.nc'];
        var_list_names       = {'fal','lai_hv','lai_lv','slhf','ssr','sp','sshf','ssrd','tp','t2m','vpd'};
        var_list_names_title = {'fal','lai\_hv','lai\_lv','slhf','ssr','sp','sshf','ssrd','tp','t2m','vpd'};
        units_list           = {'-','-','-','J.m^{-2}','J.m^{-2}','Pa','J.m^{-2}','J.m^{-2}','m','\circC','Pa'};
        lim1_var = [0.3,0.07,0.06,1e6,3e6,600,1e6,1.5e6,1e-3,5,200]; %'fal','lai_hv','lai_lv','slhf','ssr','sp','sshf','ssrd','tp','t2m','vpd'};
        
        var_list             = length(var_list_names);
        load([path_in,'GPPAnomaly_JunJul_005deg_2003_2018.mat'])
        load([path_in,'SIFAnomaly_JunJul_005deg_2003_2018.mat'])
end



for i = 1:var_list
    
    map1_flat = ncread(filename,[var_list_names{i},'_anom_r']);
    map1_flat = double(permute(map1_flat,[2,3,1]));
    % Normalise both SIF and GPP density scatter plots
    map2_flat  = GPP_anomaly_mean; % GPP anomaly
    map3_flat  = SIF_anomaly_mean; % SIF anomaly
    
    % Using the SIF data to filter Land pixels in temperature
    for a=1:size(SIF_anomaly_mean,3)
        aux = map1_flat(:,:,a);
        aux(isnan(squeeze(map2_flat(:,:,a))) | isnan(squeeze(map3_flat(:,:,a))))= NaN;
        map1_flat(:,:,a) = aux;
    end
    
    
    %------------------------------------------
    % Create the density scatter plot for SIF
    %------------------------------------------
   
    if abs(max(max(max(map1_flat))))>abs(min(min(min(map1_flat))))
            a_2     = abs(max(max(max(map1_flat))));
            n_steps = 20;            
            xedges  = [-a_2:(a_2.*2)/n_steps:a_2]; % Symmetry        
    else
            a_1     = (abs(min(min(min(map1_flat)))));
            n_steps = 20;
            xedges  = [-a_1:(a_1.*2)/n_steps:a_1]; % Symmetry
    end
    yedges_sif    = [-0.2:0.01:0.2]; % 41
    density_sif   = NaN(size(SIF_anomaly_mean,3),length(xedges)-1,length(yedges_sif)-1);
    
    
    for y = 1:size(SIF_anomaly_mean,3)
        map1 = squeeze(map1_flat(:,:,y)); map1 = map1(:);
        map3 = squeeze(map3_flat(:,:,y)); map3 = map3(:);
        % Create a 2D histogram using histogram2
        H = histogram2(map1, map3, xedges,yedges_sif); %   H = histogram2(map1, map3, numBins);
        
        title(num2str(years_covered_GPP(y)));
        % Calculate the centers of the bins
        %     x_center = (H.XBinEdges(1:end-1) + H.XBinEdges(2:end)) / 2;
        %     y_center = (H.YBinEdges(1:end-1) + H.YBinEdges(2:end)) / 2;
        % Normalize the histogram counts to get density
        density_sif(y,:,:) = H.Values / sum(H.Values(:));
        [R_SIF(i,y,:,:),P_SIF(i,y,:,:)] = corrcoef(map1,map3,'rows','complete');
    end
    
    
    f = figure;
    f.Position = [100 100 3000 2000];
    for y = 1:size(SIF_anomaly_mean,3)
        % Create the density scatter plot
        subplot(4,4,y)
        h = imagesc(xedges, yedges_sif,squeeze(density_sif(y,:,:)));      %imagesc(x_center(:), y_center(:),squeeze(density_sif(y,:,:)));
        set(gca, 'YDir', 'normal');
        hold on
        plot(xedges,zeros(1,length(xedges)),'w--')
        plot(zeros(1,length(yedges_sif)),yedges_sif,'w--')
        colormap('jet');
        h = colorbar;
        ylabel(h,'normalized density')
        ylabel('SIF Anomaly');
        xlabel([var_list_names_title{i},' anomaly ']);
        set(gca,'Fontsize',18)
        title(num2str(years_covered_GPP(y)));
    end
    if save_flag
        switch years_
            case 'initial'
                saveas(gcf,[path_fig,'LH_SIF_anomaly_Arctic_005deg_2002_2016.png']);
                saveas(gcf,[path_fig,'LH_SIF_anomaly_Arctic_005deg_2002_2016.fig'])
                
                
            case 'final'
                saveas(gcf,[path_fig,'LH_SIF_anomaly_Arctic_005deg_2003_2018.png']);
                saveas(gcf,[path_fig,'LH_SIF_anomaly_Arctic_005deg_2003_2018.fig'])
        end
    end
    
    
    %------------------------------------------
    % Density scatter plots for GPP
    %------------------------------------------
    clear density_gpp
    
    if abs(max(max(max(map1_flat))))>abs(min(min(min(map1_flat))))
            a_2     = abs(max(max(max(map1_flat))));
            n_steps = 20;            
            xedges  = [-a_2:(a_2.*2)/n_steps:a_2]; % Symmetry        
    else
            a_1     = (abs(min(min(min(map1_flat)))));
            n_steps = 20;
            xedges  = [-a_1:(a_1.*2)/n_steps:a_1]; % Symmetry
    end 
    yedges_gpp      = [-3:0.2:3]; % 31 GPP
    density_gpp     = NaN(size(GPP_anomaly_mean,3),length(xedges)-1,length(yedges_gpp)-1);
    
    for y = 1:size(SIF_anomaly_mean,3)
        map1 = squeeze(map1_flat(:,:,y)); map1 = map1(:);
        map2 = squeeze(map2_flat(:,:,y)); map2 = map2(:);
        % Create a 2D histogram using histogram2
        H = histogram2(map1, map2,xedges,yedges_gpp);%   H = histogram2(map1, map2, numBins);
        title(num2str(years_covered_GPP(y)));
        % Normalize the histogram counts to get density
        density_gpp(y,:,:) = H.Values / sum(H.Values(:));
        [R_GPP(i,y,:,:),P_GPP(i,y,:,:)] = corrcoef(map1,map2,'rows','complete');
    end
    
    % Create the density scatter plot
    f = figure;
    f.Position = [100 100 3000 2000];
    for y = 1:size(GPP_anomaly_mean,3)
        % Create the density scatter plot
        subplot(4,4,y)
        h = imagesc(xedges, yedges_gpp,squeeze(density_gpp(y,:,:))); %     imagesc(x_center(:), y_center(:),squeeze(density_gpp(y,:,:)));
        set(gca, 'YDir', 'normal');
        hold on
        plot(xedges,zeros(1,length(xedges)),'w--')
        plot(zeros(1,length(yedges_gpp)),yedges_gpp,'w--')
        colormap('jet');
        h = colorbar;
        ylabel(h,'normalized density')
        ylabel('GPP Anomaly');
        xlabel([var_list_names_title{i},' anomaly ']);
        set(gca,'Fontsize',18)
        title(num2str(years_covered_GPP(y)));
    end
    
    if save_flag
        switch years_
            case 'initial'
                saveas(gcf,[path_fig,'LH_GPP_anomaly_Arctic_005deg_2002_2016.png']);
                saveas(gcf,[path_fig,'LH_GPP_anomaly_Arctic_005deg_2002_2016.fig'])
                
                
            case 'final'
                saveas(gcf,[path_fig,'LH_GPP_anomaly_Arctic_005deg_2003_2018.png']);
                saveas(gcf,[path_fig,'LH_GPP_anomaly_Arctic_005deg_2003_2018.fig'])
        end
    end
    
    
    
    
    %------------------------------------------
    % 3-D histograms of boths SIF and GPP
    %------------------------------------------
    
    f = figure;
    f.Position = [100 100 3000 2000];
    
    yedges_gpp_norm = [-1:0.1:1];
    yedges_sif_norm = [-1:0.1:1];
    
    for y = 1:size(SIF_anomaly_mean,3)
        map1 = squeeze(map1_flat(:,:,y)); map1 = map1(:);
        map2 = squeeze(map2_flat(:,:,y)); map2 = map2(:);
        
        % In order to normalize so the 0 still is the 0 in the original series;
        % both min and max shold be simetric
        min_map2 = -4; max_map2 = 4;
        min_map3 = -0.2; max_map3 =0.2;
        
        % Apply normalization formula
        norm_map2 = ((map2 - min_map2) ./ (max_map2 - min_map2)).* 2 - 1;
        map3 = squeeze(map3_flat(:,:,y)); map3 = map3(:);
        norm_map3 = ((map3 - min_map3) ./ (max_map3 - min_map3)).* 2 - 1;
        
        
        %     norm_map2 = rescale(map2,-1,1,'InputMin',-3,'InputMax',3);
        %     norm_map3 = rescale(map3,-1,1,'InputMin',-0.2,'InputMax',0.2);
        
        subplot(4,4,y)
        % Create a 2D histogram using histogram2
        %   H1 = histogram2(map1, norm_map2, numBins,'Normalization','probability','DisplayStyle','bar3','FaceAlpha',0.6);
        H1 = histogram2(map1, norm_map2, xedges,yedges_gpp_norm,'Normalization','probability','DisplayStyle','bar3','FaceAlpha',0.6);
        
        
        H1 = H1;
        hold on
        
        %     H2 = histogram2(map1, norm_map3, numBins,'Normalization','probability','DisplayStyle','bar3','FaceAlpha',0.3);
        H2 = histogram2(map1, norm_map3,xedges, yedges_sif_norm,'Normalization','probability','DisplayStyle','bar3','FaceAlpha',0.3);
        
        H2 = H2;
        set(gca,'fontsize',16)
        legend('GPP anomaly','SIF anomaly','Location','best')
        xlabel([var_list_names_title{i},' anomaly  [',units_list{i},']' ]);
        ylabel('Norm. anomaly')
        zlabel('Norm. frequency')
        title(num2str(years_covered_GPP(y)));
        
        % Normalize the histogram counts to get density
        density_gpp_norm(y,:,:) = H1.Values / sum(H1.Values(:));
        density_sif_norm(y,:,:) = H2.Values / sum(H2.Values(:));
        
        % compute the correlation between sif and the variable and gpp and
        % the variable
        % Compute Pearson correlation coefficient
%         [corr_coef_sif(i,y), p_val_sif(i,y)] = custom_corr(map1,norm_map3);
        [R_sif(i,y,:,:),P_sif(i,y,:,:)] = corrcoef(map1,norm_map3,'rows','complete');
        
        % Compute Pearson correlation coefficient
%         [corr_coef_gpp(i,y), p_val_gpp(i,y)] = custom_corr(map1,norm_map2);
        [R_gpp(i,y,:,:),P_gpp(i,y,:,:)] = corrcoef(map1,norm_map2,'rows','complete');

        
        
    end
    
    if save_flag
        switch years_
            case 'initial'
                saveas(gcf,[path_fig,'Hist_LH_SIFGPP_norm_anomaly_Arctic_005deg_2002_2016.png']);
                saveas(gcf,[path_fig,'Hist_LH_GPP_anomaly_Arctic_005deg_2002_2016.fig'])
                
                
            case 'final'
                saveas(gcf,[path_fig,'Hist_LH_GPP_anomaly_Arctic_005deg_2003_2018.png']);
                saveas(gcf,[path_fig,'Hist_LH_GPP_anomaly_Arctic_005deg_2003_2018.fig'])
        end
    end
    


    if plotting
        %------------------------------------------
        % Normalise both SIF and GPP density scatter plots in 2D
        % Create the density scatter plot
        %------------------------------------------
        f = figure;
        f.Position = [100 100 3000 2000];
        
        x_center = (xedges(1:end-1) + xedges(2:end)) / 2;
        y_center = (yedges_gpp_norm(1:end-1) + yedges_gpp_norm(2:end)) / 2;    % are normalized, i.e. same for SIF
        
        numRows = 4;
        numCols = 4;
        subplotWidth = 0.13; % width of each subplot
        subplotHeight = 0.15; % height of each subplot
        xSpacing = 0.09; % horizontal spacing between subplots
        ySpacing = 0.07; % vertical spacing between subplots
        
        for y = 1:size(GPP_anomaly_mean,3)
            row = floor((y-1) / numCols);
            col = mod((y-1), numCols);
            left = (col * (subplotWidth + xSpacing))+xSpacing;
            bottom = 1 - (row + 1) * (subplotHeight + ySpacing);
            
            ax1 = subplot('Position', [left, bottom, subplotWidth, subplotHeight]);
            %             ax1 =subplot(4,4,y);
            contourf(ax1,x_center,y_center,squeeze(density_gpp_norm(y,:,:))',10);% 20 levels
            colormap(ax1, flipud(bone));  % Set 'bone' colormap for the first axes
            set(ax1, 'FontSize', 18);  % Set font size for the first axes
            hold on;
            % Get the position of the subplot
            pos = get(ax1, 'Position');
            ax2 = axes('Position', pos);
            contour(ax2,x_center,y_center,squeeze(density_sif_norm(y,:,:))',10, 'LineWidth', 3);
            colormap(ax2, jet);  % Set 'jet' colormap for the second axes
            ax2.Color = 'none';  % Make the second axes transparent
            set(ax2, 'FontSize', 18);  % Set font size for the second axes
            hold on
            plot(x_center(:),zeros(1,length(x_center(:))),'k--')
            set(gca,'Fontsize',18)
            plot(zeros(1,length(y_center(:))),y_center(:),'k--')
            set(gca,'Fontsize',18)
            xlabel(ax2,[var_list_names_title{i},' anomaly [',units_list{i},']' ]);
            title(num2str(years_covered_GPP(y)));
            
            % Link the axes
            linkaxes([ax1, ax2]);
            
            % Set the second axes to be on top
            ax2.Position = ax1.Position;
            
            
            % Adjust the Y-axis limit of the first axes
            ylim(ax1, [-0.4, 0.4]);
            ylim(ax2, [-0.4, 0.4]);
            
            % Set the ylabel for the first axes
            ylabel(ax1, {'GPP & SIF'; 'norm. Anomaly'});
            
            xlim([-lim1_var(i),lim1_var(i)])
            
            if sum(squeeze(P_gpp(i,y,1,2)) && squeeze(P_sif(i,y,1,2)))==0
                % Display the correlation coefficient on the plot
                annotation_text = sprintf('r_{GPP} = %.2f            r_{SIF} = %.2f', squeeze(R_gpp(i,y,1,2)),squeeze(R_sif(i,y,1,2)));
                text(ax1, 'Units', 'normalized', 'Position', [0.05, 0.08], 'String', annotation_text, 'FontSize', 18, 'Color', 'k', 'FontWeight', 'bold','BackgroundColor', [1, 1, 1, 0.9], ... % white background with 50% opacity
           'EdgeColor', 'none');
            else
                
                % Display the correlation coefficient on the plot
                annotation_text = sprintf('r_{GPP} = %.2f p_{value-GPP} <0.01 \nr_{SIF} = %.2f p_{value-SIF} = <0.01', squeeze(R_gpp(i,y,1,2)),squeeze(P_gpp(i,y,1,2)),squeeze(R_sif(i,y,1,2)),squeeze(P_sif(i,y,1,2)));
                text(ax1, 'Units', 'normalized', 'Position', [0.05, 0.08], 'String', annotation_text, 'FontSize', 18, 'Color', 'k', 'FontWeight', 'bold');
            end
            
            
        end
        
        
        % Add colorbars after the axes are correctly positioned
        cb1 = colorbar(ax1, 'Position', [pos(1)+pos(3)+0.01, pos(2), 0.02, pos(4)], 'FontSize', 18,'FontWeight','bold');
        ylabel(cb1,{'GPP norm. density'})
        cb2 = colorbar(ax2, 'southoutside', 'FontSize', 18,'FontWeight','bold');
        ylabel(cb2,{'SIF norm. density'})
        % Adjust the position of the second colorbar
        set(cb2, 'Position', [pos(1), pos(2)-0.065, pos(3), 0.02]);
        
        
        
%         if save_flag
            switch years_
                case 'initial'
                    
                case 'final'
                    saveas(gcf,[path_fig,filesep,var_list_names{i},'_SIFGPP_norm_anomaly_Arctic_005deg_2003_2018.png']);
                    saveas(gcf,[path_fig,filesep,var_list_names{i},'_SIFGPP_norm_anomaly_Arctic_005deg_2003_2018.fig']);
                    print([path_fig,filesep,var_list_names{i},'_SIFGPP_norm_anomaly_Arctic_005deg_2003_2018.eps'], '-depsc', '-r600', '-opengl');
            end
%         end
        
        
        %------------------------------------------
        % SIF and GPP density scatter plots in 2D
        % Create the density scatter plot different y axis
        %------------------------------------------
        
        f = figure;
        f.Position = [100 100 2000 2400];
        
        x_center = (xedges(1:end-1) + xedges(2:end)) / 2;
        y_center_gpp = (yedges_gpp(1:end-1) + yedges_gpp(2:end)) / 2;
        y_center_sif = (yedges_sif(1:end-1) + yedges_sif(2:end)) / 2;
        
        numRows = 4;
        numCols = 4;
        subplotWidth = 0.13; % width of each subplot
        subplotHeight = 0.15; % height of each subplot
        xSpacing = 0.09; % horizontal spacing between subplots
        ySpacing = 0.07; % vertical spacing between subplots
        
        for y = 1:size(GPP_anomaly_mean, 3)
            row = floor((y-1) / numCols);
            col = mod((y-1), numCols);
            left = (col * (subplotWidth + xSpacing))+xSpacing;
            bottom = 1 - (row + 1) * (subplotHeight + ySpacing);
            
            ax1 = subplot('Position', [left, bottom, subplotWidth, subplotHeight]);
            contourf(ax1, x_center, y_center_gpp, squeeze(density_gpp(y,:,:))', 10, 'LineStyle', 'none'); % 20 levels
            colormap(ax1, flipud(bone)); % Set 'bone' colormap for the first axes
            set(ax1, 'FontSize', 18); % Set font size for the first axes
            hold on;
            
            % Overlay the second axes on top of the first
            ax2 = axes('Position', ax1.Position);
            contour(ax2, x_center, y_center_sif, squeeze(density_sif(y,:,:))', 10, 'LineWidth', 1);
            colormap(ax2, jet); % Set 'jet' colormap for the second axes
            set(ax2, 'FontSize', 18); % Set font size for the second axes
            hold on;
            
            % Ensure the second axes is transparent and positioned correctly
            ax2.Color = 'none';
            ax2.XColor = 'none';
            ax2.YAxisLocation = 'right';
            ax2.Position = ax1.Position;
            title(ax2, num2str(years_covered_GPP(y)));
            
            % Plot the zero lines on both axes
            plot(ax1, x_center(:), zeros(1, length(x_center(:))), 'k--')
            plot(ax1, zeros(1, length(y_center_gpp(:))), y_center_gpp(:), 'k--')
            % plot(ax2, x_center(:), zeros(1, length(x_center(:))), 'k--')
            % plot(ax2, zeros(1, length(y_center_sif(:))), y_center_sif(:), 'k--')
            
            % Set labels and titles
            xlabel(ax1, [var_list_names_title{i}, ' anomaly [', units_list{i}, ']']);
            ylabel(ax1, 'GPP anomaly');
            ylabel(ax2, 'SIF anomaly');
            
            % Set x and y limits
            xlim(ax1, [-lim1_var(i), lim1_var(i)])
            xlim(ax2, [-lim1_var(i), lim1_var(i)])
            % ylim(ax1, [y_center_gpp(1)+0.2, y_center_gpp(end)-0.2])
            % ylim(ax2, [y_center_sif(1), y_center_sif(end)])
            ylim(ax1, [-2.0, 2.0])
            ylim(ax2, [-0.1, 0.1])
            
            
            if sum(squeeze(P_gpp(i,y,1,2)) && squeeze(P_sif(i,y,1,2)))==0
                % Display the correlation coefficient on the plot
                annotation_text = sprintf('r_{GPP} = %.2f         r_{SIF} = %.2f', squeeze(R_GPP(i,y,1,2)),squeeze(R_SIF(i,y,1,2)));
                text(ax1, 'Units', 'normalized', 'Position', [0.05, 0.08], 'String', annotation_text, 'FontSize', 18, 'Color', 'k', 'FontWeight', 'bold');
            else
                
                % Display the correlation coefficient on the plot
                annotation_text = sprintf('r_{GPP} = %.2f p_{value-GPP} <0.01 \nr_{SIF} = %.2f p_{value-SIF} = <0.01', squeeze(R_GPP(i,y,1,2)),squeeze(P_GPP(i,y,1,2)),squeeze(R_SIF(i,y,1,2)),squeeze(P_SIF(i,y,1,2)));
                text(ax1, 'Units', 'normalized', 'Position', [0.05, 0.08], 'String', annotation_text, 'FontSize', 18, 'Color', 'k', 'FontWeight', 'bold');
            end
        end
        
        % Get the position of the subplot
        pos = get(ax1, 'Position');
        % Add colorbars after the axes are correctly positioned
        cb1 = colorbar(ax1, 'Position', [pos(1)+pos(3)+0.06, pos(2), 0.02, pos(4)], 'FontSize', 18,'FontWeight','bold');
        ylabel(cb1,{'GPP density'})
        cb2 = colorbar(ax2, 'southoutside', 'FontSize', 18,'FontWeight','bold');
        ylabel(cb2,{'SIF density'})
        % Adjust the position of the second colorbar
        set(cb2, 'Position', [pos(1), pos(2)-0.065, pos(3), 0.02]);
        
        
        
%         if save_flag
            switch years_
                case 'initial'
                    
                    
                case 'final'
                    saveas(gcf,[path_fig,filesep,var_list_names{i},'_SIFGPP_anomaly_Arctic_005deg_2003_2018.png']);
                    saveas(gcf,[path_fig,filesep,var_list_names{i},'_SIFGPP_anomaly_Arctic_005deg_2003_2018.fig']);
                    print([path_fig,filesep,var_list_names{i},'_SIFGPP_anomaly_Arctic_005deg_2003_2018.eps'], '-depsc', '-r600', '-opengl');
            end
%         end
        
    end
    close all
    
end % varlist



%% ------------------------------------------
%  Plotting the correlation coefficients
%------------------------------------------
correlationCoefficients = squeeze(R_GPP(4:end,:,1,2));
years = 2003:2018; % Replace with your actual years
variableNames =  {'slhf','ssr','sp','sshf','ssrd','tp','t2m','vpd'};

f = figure;
f.Position = [100 100 800 500];
        
% Create the heatmap
cmap = redblue;
h = heatmap(correlationCoefficients, 'XLabel', 'Year (Jun-Jul)',...
            'Title', 'Correlation Coefficients - r(var anom,GPP anom)', 'Colormap', cmap);

% Set the custom x and y axis labels
h.XDisplayLabels = string(years);
h.YDisplayLabels = variableNames;
h.CellLabelFormat = '%.2f';
set(gca,'Fontsize',18)
% caxis([-0.8,0.8])
caxis([-1,1])

% if save_flag
    switch years_
        case 'initial'
            
            
        case 'final'
            saveas(gcf,[path_fig,filesep,'Corr_GPP_heatmap_2003_2018.png']);
            saveas(gcf,[path_fig,filesep,'Corr_GPP_heatmap_2003_2018.fig']);
            print([path_fig,filesep,'Corr_GPP_heatmap_2003_2018.eps'], '-depsc', '-r600', '-opengl');
    end
% end

correlationCoefficients = squeeze(R_SIF(4:end,:,1,2));
years = 2003:2018; % Replace with your actual years
variableNames =  {'slhf','ssr','sp','sshf','ssrd','tp','t2m','vpd'};

% Create the heatmap
f = figure;
f.Position = [100 100 800 500];
        
h = heatmap(correlationCoefficients, 'XLabel', 'Year (Jun-Jul)',...
            'Title', 'Correlation Coefficients - r(var anom,SIF anom) ', 'Colormap', cmap);

% Set the custom x and y axis labels
h.XDisplayLabels = string(years);
h.YDisplayLabels = variableNames;
h.CellLabelFormat = '%.2f';
set(gca,'Fontsize',18)
caxis([-1,1])

% if save_flag
    switch years_
        case 'initial'
            
            
        case 'final'
            saveas(gcf,[path_fig,filesep,'Corr_SIF_heatmap_2003_2018.png']);
            saveas(gcf,[path_fig,filesep,'Corr_SIF_heatmap_2003_2018.fig']);
            print([path_fig,filesep,'Corr_SIF_heatmap_2003_2018.eps'], '-depsc', '-r600', '-opengl');
    end
% end

%% ------------------------------------------
% Normalise both SIF and GPP density scatter plots in 2D
% Create the density scatter plot
%------------------------------------------
% Get the normalised density scatter plots in 1 year = 2016

clear density_gpp_norm density_sif_norm


for i = 4:var_list
    map1_flat = ncread(filename,[var_list_names{i},'_anom_r']);
    map1_flat = double(permute(map1_flat,[2,3,1]));
    % Normalise both SIF and GPP density scatter plots
    map2_flat  = GPP_anomaly_mean; % GPP anomaly
    map3_flat  = SIF_anomaly_mean; % SIF anomaly
    
    % Using the SIF data to filter Land pixels in temperature
    for a=1:size(SIF_anomaly_mean,3)
        aux = map1_flat(:,:,a);
        aux(isnan(squeeze(map2_flat(:,:,a))) | isnan(squeeze(map3_flat(:,:,a))))= NaN;
        map1_flat(:,:,a) = aux;
    end
    
    
    %------------------------------------------
    % Create the density scatter plot for SIF
    %------------------------------------------
    %     xedges        = [-5:0.5:5]; % 21
    if abs(max(max(max(map1_flat))))>abs(min(min(min(map1_flat))))
            a_2     = abs(max(max(max(map1_flat))));
            n_steps = 20;            
            xedges  = [-a_2:(a_2.*2)/n_steps:a_2]; % Symmetry        
    else
            a_1     = (abs(min(min(min(map1_flat)))));
            n_steps = 20;
            xedges  = [-a_1:(a_1.*2)/n_steps:a_1]; % Symmetry
    end
    
    
    %------------------------------------------
    % 3-D histograms of boths SIF and GPP
    %------------------------------------------
   
    
    yedges_gpp_norm = [-1:0.1:1];
    yedges_sif_norm = [-1:0.1:1];
    
    
    y = 14 % 2016
    map1 = squeeze(map1_flat(:,:,y)); map1 = map1(:);
    map2 = squeeze(map2_flat(:,:,y)); map2 = map2(:);
    
    % In order to normalize so the 0 still is the 0 in the original series;
    % both min and max shold be simetric
    min_map2 = -4; max_map2 = 4;
    min_map3 = -0.2; max_map3 =0.2;
    
    % Apply normalization formula
    norm_map2 = ((map2 - min_map2) ./ (max_map2 - min_map2)).* 2 - 1;
    map3 = squeeze(map3_flat(:,:,y)); map3 = map3(:);
    norm_map3 = ((map3 - min_map3) ./ (max_map3 - min_map3)).* 2 - 1;
    
    % Create a 2D histogram using histogram2
    figure
    H1 = histogram2(map1, norm_map2, xedges,yedges_gpp_norm,'Normalization','probability','DisplayStyle','bar3','FaceAlpha',0.6);
    H1 = H1;
    density_gpp_norm(1,:,:) = H1.Values / sum(H1.Values(:));     % Normalize the histogram counts to get density
    close(gcf)
    
    figure
    H2 = histogram2(map1, norm_map3,xedges, yedges_sif_norm,'Normalization','probability','DisplayStyle','bar3','FaceAlpha',0.3);
    H2 = H2; 
    density_sif_norm(1,:,:) = H2.Values / sum(H2.Values(:));
    close(gcf)    
    
     % compute the correlation between sif and the variable and gpp and
     % the variable
     % Compute Pearson correlation coefficient
%      [corr_coef_sif_year(i), p_val_sif_year(i)] = custom_corr(map1,norm_map3);
     [R_sif(i,:,:),P_sif(i,:,:)] = corrcoef(map1,norm_map3,'rows','complete');

     
     % Compute Pearson correlation coefficient
%      [corr_coef_gpp_year(i), p_val_gpp_year(i)] = custom_corr(map1,norm_map2);
     [R_gpp(i,:,:),P_gpp(i,:,:)] = corrcoef(map1,norm_map2,'rows','complete');

 
    
    x_center = (xedges(1:end-1) + xedges(2:end)) / 2;
    y_center = (yedges_gpp_norm(1:end-1) + yedges_gpp_norm(2:end)) / 2;    % are normalized, i.e. same for SIF
    
    f = figure(200);
    f.Position = [100 100 3000 2000];
    ax1 =subplot(2,4,i-3);
    contourf(ax1,x_center,y_center,squeeze(density_gpp_norm(1,:,:))',10);% 20 levels
    colormap(ax1, flipud(bone));  % Set 'bone' colormap for the first axes
    set(ax1, 'FontSize', 18);  % Set font size for the first axes
    hold on;
    % Get the position of the subplot
    pos = get(ax1, 'Position');
    ax2 = axes('Position', pos);
    contour(ax2,x_center,y_center,squeeze(density_sif_norm(1,:,:))',10, 'LineWidth', 3);
    colormap(ax2, jet);  % Set 'jet' colormap for the second axes
    ax2.Color = 'none';  % Make the second axes transparent
    set(ax2, 'FontSize', 18);  % Set font size for the second axes
    hold on
    plot(x_center(:),zeros(1,length(x_center(:))),'k--')
    set(gca,'Fontsize',18)
    plot(zeros(1,length(y_center(:))),y_center(:),'k--')
    set(gca,'Fontsize',18)
    xlabel(ax2,[var_list_names_title{i},' anomaly [',units_list{i},']' ]);
    title([var_list_names_title{i},' - anomalies ',num2str(years_covered_GPP(y))]);
    
    % Link the axes
    linkaxes([ax1, ax2]);
    
    % Set the second axes to be on top
    ax2.Position = ax1.Position;
    
    
    % Adjust the Y-axis limit of the first axes
    ylim(ax1, [-0.3, 0.4]);
    ylim(ax2, [-0.3, 0.4]);
    
    % Set the ylabel for the first axes
    ylabel(ax1, {'GPP & SIF'; 'norm. Anomaly'});
    
    xlim([-lim1_var(i),lim1_var(i)])
    
    if sum(squeeze(P_gpp(i,1,2)) && squeeze(P_sif(i,1,2)))==0
        % Display the correlation coefficient on the plot
        annotation_text = sprintf('r_{GPP} = %.2f p_{value-GPP} <0.01 \nr_{SIF} = %.2f p_{value-SIF} <0.01', squeeze(R_gpp(i,1,2)),squeeze(R_sif(i,1,2)));
        text(ax1, 'Units', 'normalized', 'Position', [0.05, 0.08], 'String', annotation_text, 'FontSize', 18, 'Color', 'k', 'FontWeight', 'bold');
    else
        
        % Display the correlation coefficient on the plot
        annotation_text = sprintf('r_{GPP} = %.2f p_{value-GPP} <0.01 \nr_{SIF} = %.2f p_{value-SIF} = <0.01', squeeze(R_gpp(i,1,2)),squeeze(P_gpp(i,1,2)),squeeze(R_sif(i,1,2)),squeeze(P_sif(i,1,2)));
        text(ax1, 'Units', 'normalized', 'Position', [0.05, 0.08], 'String', annotation_text, 'FontSize', 18, 'Color', 'k', 'FontWeight', 'bold');
    end
    
    
end   

% Add colorbars after the axes are correctly positioned
cb1 = colorbar(ax1, 'Position', [pos(1)+pos(3)+0.01, pos(2), 0.02, pos(4)], 'FontSize', 18);
ylabel(cb1,{'GPP norm. density'})
cb2 = colorbar(ax2, 'southoutside', 'FontSize', 18);
ylabel(cb2,{'SIF norm. density'})
% Adjust the position of the second colorbar
set(cb2, 'Position', [pos(1), pos(2)-0.065, pos(3), 0.02]);


% 
% if save_flag
    switch years_
        case 'initial'
            
        case 'final'
            saveas(gcf,[path_fig,filesep,'climate_SIFGPP_norm_anomaly_Arctic_005deg_2016.png']);
            saveas(gcf,[path_fig,filesep,'climate_SIFGPP_norm_anomaly_Arctic_005deg_2016.fig']);
            print([path_fig,filesep,'climate_SIFGPP_norm_anomaly_Arctic_005deg_2016.eps'], '-depsc', '-r600', '-opengl');
    end
% end
















% #######################################################################################
% #######################################################################################
% #######################################################################################
% #######################################################################################
% #######################################################################################
