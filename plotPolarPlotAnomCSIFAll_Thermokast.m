% ###################################################################
% ################# plotPolarPlotAnomCSIFAll_Thermokast   ###########
% This script reads the cSIF database anomalies from 2003-2016 (or 2003-2018) and computes 
% the OLS regression with the corresponding airT anomalies for the same
% years. 
% We also use the Thermoskast masps available in 
% https://daac.ornl.gov/cgi-bin/dsviewer.pl?ds_id=1332 and
% published in https://www.nature.com/articles/ncomms13043
% to classify the slopes according to the Thermokast presence/degree.
%
% Data was extracted from the cSIF database generated with OCO-2 SIF data and 
% MODIS surface reflectance every 4-day at 0.05 degrees. 
% This database also disentangles the sunlit and all sky SIF prediction.
% 
% %Zhang, Y., Joiner, J., Alemohammad, S. H., Zhou, S., & Gentine, P. (2018).
% A global spatially contiguous solar-induced fluorescence (CSIF) dataset 
% using neural networks. Biogeosciences, 15(19), 5779-5800.
% https://bg.copernicus.org/articles/15/5779/2018/

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
% airTemperatureswitch years_
switch years_
    case 'cSIF_all'
        load([path_in,'airTAnomaly_JunJul_005deg_2003_2016.mat']) % 'years_cov','t2m_summer_anomaly_year','lon','lat'))
        t2m_anomaly_mean = permute(t2m_anomaly_mean,[2,3,1]);     % change the dimensions long lat years
        
    case 'cSIF_clear'
        load([path_in,'airTAnomaly_JunJul_005deg_2003_2018.mat']) % 'years_cov','t2m_summer_anomaly_year','lon','lat')) 
end

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



% GPP 
load([path_in,'SIF_GPP_Z_score_TrendArticSIFGPPFluxSat_JJ_2003_2018.mat'],'GPP_TrendStore')  
switch years_
    case 'cSIF_all'
        GPP_TrendStore     = GPP_TrendStore(:,:,1:14);% trim the dataseries up to 2016
end
GPP_mean           = mean(GPP_TrendStore,3,'omitnan');
GPP_median         = median(GPP_TrendStore,3,'omitnan');
GPP_anomaly_mean   = GPP_TrendStore-GPP_mean;
GPP_anomaly_median = GPP_TrendStore-GPP_median;



switch years_
    
    case 'cSIF_all'
        % cSIF data
        load([path_in,'cSIFAnomaly_JunJul_05deg_2003_2016.mat'],'cSIF_anomaly_mean','cSIF_anomaly_median',...
            'cSIF_median','cSIF_mean','longitude','latitude','years_covered_cSIF')
        cSIF_anomaly_mean  = permute(cSIF_anomaly_mean,[2,3,1]);   % change the dimensions long lat years
        
        
    case 'cSIF_clear'
        load([path_in,'cSIFclearAnomaly_JunJul_005deg_2003_2018.mat'],'cSIF_anomaly_mean','cSIF_anomaly_median',...
            'cSIF_median','cSIF_mean','longitude','latitude','years_covered_cSIF')
        cSIF_anomaly_mean  = permute(cSIF_anomaly_mean,[2,3,1]);   % change the dimensions long lat years
end
%% --------------------------
% Plotting and computing the slope
% ---------------------------


% if plotting == 1
    f = figure;
    f.Position = [100 100 2000 1200];
    
    n = 5; % Number of colors according to the Thermokast classification
    % Define custom colormap with high contrast
    colors = [227, 227, 227;182, 146, 222;131, 89, 179;110, 46, 184;39, 4, 79]./255;
    tickLabels = {'None','Low','Moderate','High','Very High'}; 
        
    % GPP anomalies scatter plots vs. Air T anomalies
    f = figure;
    f.Position = [100 100 2000 1200];
    for j=1:5 % Number of classes
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
        subplot(2,3,j)
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
        x_new = [-5:0.1:5]';
        y = auxGPP_an_store(~isnan(auxGPP_an_store) & ~isnan(auxT_an_store));
        X = [ones(length(x),1) x];
        
        % matlab OLS linear regressor
        b          = X\y;
        GPP_B(j,:) = b;
        yCalc2     = X*b;       
        Rsq2       = 1 - sum((y - yCalc2).^2)/sum((y - mean(y)).^2);
        
        % Compute other regressors info
        mdl                 = fitlm(x,y); % OLS linear
        slopeGPP{j}         = mdl.Coefficients.Estimate(2); % Slope (coefficient of X)
        slopeGPP_SE{j}      = mdl.Coefficients.SE(2);       % Standard error of the slope
        resultsGPP{j}       = mdl;
        slope_conf_interval = coefCI(mdl);
        slopeGPP_CI{j,:}    = slope_conf_interval(2,:); 
        P(j)  = mdl.Coefficients.pValue(2);
        fprintf('pValue = %.2e\n', P(j));
        
     
        % Plot the fitting with the confidence levels 
        % plot(x,yCalc2,'k')
        [GPP_pred, GPP_pred_CI] = predict(mdl, x_new);
        plot(x_new, GPP_pred, '-r', 'LineWidth', 2);
        % Plot confidence intervals as shaded areas
        fill([x_new; flipud(x_new)], [GPP_pred_CI(:,1); flipud(GPP_pred_CI(:,2))], ...
         'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none');  % Shaded area for confidence intervals

        
        
        % Matlab linear regressor SVM
        [Mdl,FitInfo] = fitrlinear(x,y);  % svm learner        
        GPP_Bsvm(j,:) = Mdl.Beta;
        
        
        % simpleRegression Toolbox for linear regressors
        % ------------------------------------------------
        rate = 1; %[0.05 0.1 0.2 0.3 0.4 0.5 0.6]
        % Fix seed random generator (important: disable when doing the 100 realizations loop!)
        % rand('seed',12345);
        % randn('seed',12345);
        % rng(0);
        [n,~] = size(x);                 % samples x bands
        r = randperm(n);                 % random index
        ntrain = round(rate*n);          % #training samples
        Xtrain = x(r(1:ntrain),:);       % training set
        Ytrain = y(r(1:ntrain),:);       % observed training variable
        Xtest  = x(r(ntrain+1:end),:);   % test set
        Ytest  = y(r(ntrain+1:end),:);   % observed test variable
        [ntest do] = size(Ytest);
        
        METHODS = {'RLR' 'LASSO' 'ENET'}; % LINEAR  
        numModels = numel(METHODS);
        
        for m=1:numModels
            fprintf(['Training ' METHODS{m} '... \n'])
            t=cputime;
            eval(['model = train' METHODS{m} '(Xtrain,Ytrain);']); % Train the model
%             eval(['Yp = test' METHODS{m} '(model,Xtest);']);       % Test the model
%             RESULTS(m)  = assessment(Ytest, Yp, 'regress');  % assessregres(Ytest,Yp);
%             CPUTIMES(m) = cputime - t;
            MODELS{m} = model;
%             YPREDS(:,m) = Yp;
        end
        
        GPP_Brlr(j,:)   = MODELS{1}.W(2);
        GPP_Blasso(j,:) = MODELS{2}.B(:,MODELS{2}.S.Index1SE);
        GPP_Benet(j,:)  = MODELS{3}.B(:,MODELS{3}.S.Index1SE);
        clear MODELS
        % -------------------------------------------------
        
        
        
        box on
        grid on
        grid minor
        set(gca,'TickLength',[0.05, 0.01])
        set(gca,'XMinorTick','on','YMinorTick','on')
        text(-4,3,['R^2 =',num2str(round(Rsq2,2))],'Fontsize',18)
        text(-4,4,['y = ',num2str(round(b(1),6)),' + ',num2str(round(b(2),4)),'x'],'Fontsize',18)
    end
    
    if save_flag
        switch years_
            case 'initial'
                saveas(gcf,[path_fig,'Thermokast_GPP_class_scatter_2002_2016.png'])
                saveas(gcf,[path_fig,'Thermokast_GPP_class_scatter_2002_2016.fig'])
                
            case 'final'
                saveas(gcf,[path_fig,'Thermokast_GPP_class_scatter_2003_2018.png'])
                saveas(gcf,[path_fig,'Thermokast_GPP_class_scatter_2003_2018.fig'])
                
            case 'cSIF_all'
                saveas(gcf,[path_fig,'Thermokast_GPP_class_scatter_2003_2016.png'])
                saveas(gcf,[path_fig,'Thermokast_GPP_class_scatter_2003_2016.fig'])
                
            case 'cSIF_clear'
                saveas(gcf,[path_fig,'Thermokast_GPPclear_class_scatter_2003_2018.png'])
                saveas(gcf,[path_fig,'Thermokast_GPPclear_class_scatter_2003_2018.fig'])
                
        end
    end
    
    
%%   cSIF anomalies 
    f = figure;
    f.Position = [100 100 2000 1200];
    for j=1:5
        auxT_an_store   = [];
        auxSIF_an_store = [];
        
        for i=1:size(GPP_TrendStore,3)
            auxSIF = squeeze(cSIF_anomaly_mean(:,:,i));
            auxT   = squeeze(t2m_anomaly_mean(:,:,i));
            
            auxT   = auxT(CLASS_int==j);
            auxSIF = auxSIF(CLASS_int==j);
            
            auxSIF_an_store = [auxSIF_an_store(:);auxSIF(:)]; % all years
            auxT_an_store   = [auxT_an_store(:);auxT(:)]; % all years
        end
        subplot(2,3,j)
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
        Rsq2 = 1 - sum((y - yCalc2).^2)/sum((y - mean(y)).^2);
        
        % Compute other regressors info
        mdl   = fitlm(x,y); % OLS linear
        
        slopeSIF{j}         = mdl.Coefficients.Estimate(2); % Slope (coefficient of X)
        slopeSIF_SE{j}      = mdl.Coefficients.SE(2);       % Standard error of the slope
        resultsSIF{j}       = mdl;
        slope_conf_interval = coefCI(mdl);
        slopeSIF_CI{j,:}    = slope_conf_interval(2,:);  
        P(j)  = mdl.Coefficients.pValue(2);
        fprintf('pValue = %.2e\n', P(j));
        
        % Plot the fitting with the confidence levels 
        % plot(x,yCalc2,'k')
        [SIF_pred, SIF_pred_CI] = predict(mdl, x_new);
        plot(x_new, SIF_pred, '-r', 'LineWidth', 2);
        % Plot confidence intervals as shaded areas
        fill([x_new; flipud(x_new)], [SIF_pred_CI(:,1); flipud(SIF_pred_CI(:,2))], ...
         'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none');  % Shaded area for confidence intervals

        
        % Matlab linear regressor SVM
        [Mdl,FitInfo] = fitrlinear(x,y);  % svm learner        
        SIF_Bsvm(j,:) = Mdl.Beta;
        
        % simpleRegression Toolbox for linear regressors
        % ------------------------------------------------
        rate = 1; %[0.05 0.1 0.2 0.3 0.4 0.5 0.6]
        % Fix seed random generator (important: disable when doing the 100 realizations loop!)
        % rand('seed',12345);
        % randn('seed',12345);
        % rng(0);
        [n,~] = size(x);                 % samples x bands
        r = randperm(n);                 % random index
        ntrain = round(rate*n);          % #training samples
        Xtrain = x(r(1:ntrain),:);       % training set
        Ytrain = y(r(1:ntrain),:);       % observed training variable
        Xtest  = x(r(ntrain+1:end),:);   % test set
        Ytest  = y(r(ntrain+1:end),:);   % observed test variable
        [ntest do] = size(Ytest);
        
        METHODS = {'RLR' 'LASSO' 'ENET'}; % LINEAR  
        numModels = numel(METHODS);
        
        for m=1:numModels
            fprintf(['Training ' METHODS{m} '... \n'])
            t=cputime;
            eval(['model = train' METHODS{m} '(Xtrain,Ytrain);']); % Train the model
%             eval(['Yp = test' METHODS{m} '(model,Xtest);']);       % Test the model
%             RESULTS(m)  = assessment(Ytest, Yp, 'regress');  % assessregres(Ytest,Yp);
%             CPUTIMES(m) = cputime - t;
            MODELS{m} = model;
%             YPREDS(:,m) = Yp;
        end
        
        SIF_Brlr(j,:)   = MODELS{1}.W(2);
        SIF_Blasso(j,:) = MODELS{2}.B(:,MODELS{2}.S.Index1SE);
        SIF_Benet(j,:)  = MODELS{3}.B(:,MODELS{3}.S.Index1SE);
        clear MODELS
        % -------------------------------------------------
        
        box on
        grid on
        grid minor
        set(gca,'TickLength',[0.05, 0.01])
        set(gca,'XMinorTick','on','YMinorTick','on')
        text(-4,0.1,['R^2 =',num2str(round(Rsq2,2))],'Fontsize',18)
        text(-4,0.2,['y = ',num2str(round(b(1),6)),' + ',num2str(round(b(2),4)),'x'],'Fontsize',18)
    end
    
    if save_flag
        switch years_
            case 'initial'
                saveas(gcf,[path_fig,'Thermokast_cSIF_class_scatter_2002_2016.png'])
                saveas(gcf,[path_fig,'Thermokast_cSIF_class_scatter_2002_2016.fig'])
                
            case 'final'
                saveas(gcf,[path_fig,'Thermokast_cSIF_class_scatter_2003_2018.png'])
                saveas(gcf,[path_fig,'Thermokast_cSIF_class_scatter_2003_2018.fig'])
                
            case 'cSIF_all'
                saveas(gcf,[path_fig,'Thermokast_cSIF_class_scatter_2003_2016.png'])
                saveas(gcf,[path_fig,'Thermokast_cSIF_class_scatter_2003_2016.fig'])
                
            case 'cSIF_clear'
                saveas(gcf,[path_fig,'Thermokast_cSIFclear_class_scatter_2003_2018.png'])
                saveas(gcf,[path_fig,'Thermokast_cSIFclear_class_scatter_2003_2018.fig'])
        end
    end
    
%% Figure of the anomalies rates of SIF and GPP per class
% from the linear regression

% if plotting == 1
    x = [-5:0.1:5]';
    X = [ones(length(x),1) x];
    for j=1:5
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
                saveas(gcf,[path_fig,'Thermokast_SIF_linearreg_2002_2016.png'])
                saveas(gcf,[path_fig,'Thermokast_SIF_linearreg_2002_2016.fig'])
                
            case 'final'
                saveas(gcf,[path_fig,'Thermokast_cSIF_linearreg_2003_2018.png'])
                saveas(gcf,[path_fig,'Thermokast_cSIF_linearreg_2003_2018.fig'])
                
            case 'cSIF_all'
                saveas(gcf,[path_fig,'Thermokast_cSIF_linearreg_2003_2016.png'])
                saveas(gcf,[path_fig,'Thermokast_cSIF_linearreg_2003_2016.fig'])
                
            case 'cSIF_clear'
                saveas(gcf,[path_fig,'Thermokast_cSIFclear_linearreg_2003_2018.png'])
                saveas(gcf,[path_fig,'Thermokast_cSIFclear_linearreg_2003_2018.fig'])
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
                saveas(gcf,[path_fig,'Thermokast_GPP_linearreg_2002_2016.png'])
                saveas(gcf,[path_fig,'Thermokast_GPP_linearreg_2002_2016.fig'])
                
            case 'final'
                saveas(gcf,[path_fig,'Thermokast_GPP_linearreg_2003_2018.png'])
                saveas(gcf,[path_fig,'Thermokast_GPP_linearreg_2003_2018.fig'])
                
            case 'cSIF_all'
                saveas(gcf,[path_fig,'Thermokast_GPP_linearreg_2003_2016.png'])
                saveas(gcf,[path_fig,'Thermokast_GPP_linearreg_2003_2016.fig'])
                
            case 'cSIF_clear'
                saveas(gcf,[path_fig,'Thermokast_GPP_linearreg_2003_2018.png'])
                saveas(gcf,[path_fig,'Thermokast_GPP_linearreg_2003_2018.fig'])
        end
    end
    
% end



if plotting == 1
    
    figure
    for j=1:5
        yyaxis left
        h1 = plot(j,SIF_B(j,2),'o','Markersize',26,'markerfacecolor',colors(j,:),'markeredgecolor','k');
        hold on
        yyaxis right
        h2 = plot(j,GPP_B(j,2),'d','Markersize',26,'markerfacecolor',colors(j,:),'markeredgecolor','k');
        hold on
    end
    yyaxis left
    plot([1:5],SIF_B(:,2),'k-')
    yyaxis right
    plot([1:5],GPP_B(:,2),'k-')
    yyaxis right
    ylabel('GPP anomaly vs. air T anomaly')
    yyaxis left
    ylabel('SIF anomaly vs. air T anomaly')
    xticks([1:5])
    xticklabels({'None','Low','Moderate','High','Very High'}); 
    xtickangle(45)
    box on
    grid on
    grid minor
    set(gca,'Fontsize',18)
    set(gca,'TickLength',[0.05, 0.005])
    set(gca,'XMinorTick','on','YMinorTick','on')
    xlim([0,6])
    title('GPP and SIF anomalies response to air Temperature')
    lgnd = legend([h1,h2],'SIF','GPP','Location','west','orientation','horizontal');
    set(lgnd,'color',[0.8,0.8,0.8]);
    if save_flag
        switch years_
            case 'initial'
                saveas(gcf,[path_fig,'Thermokast_GPPcSIF_anom_changevsairT_2002_2016.png'])
                saveas(gcf,[path_fig,'Thermokast_GPPcSIF_anom_changevsairT_2002_2016.fig'])
                
            case 'final'
                saveas(gcf,[path_fig,'Thermokast_GPPcSIF_anom_changevsairT_2003_2018.png'])
                saveas(gcf,[path_fig,'Thermokast_GPPcSIF_anom_changevsairT_2003_2018.fig'])
            
            case 'cSIF_all'
                saveas(gcf,[path_fig,'Thermokast_GPPcSIF_anom_changevsairT_2003_2016.png'])
                saveas(gcf,[path_fig,'Thermokast_GPPcSIF_anom_changevsairT_2003_2016.fig'])
                
            case 'cSIF_clear'
                saveas(gcf,[path_fig,'Thermokast_GPPcSIFclear_anom_changevsairT_2003_2018.png'])
                saveas(gcf,[path_fig,'Thermokast_GPPcSIFclear_anom_changevsairT_2003_2018.fig'])
                
        end
    end
    
end


%% Order the data in descend order
% if plotting == 1
% using the OLS regression
%     [~,posSIF] = sort(SIF_B(:,2),'descend');
%     [~,posGPP]  = sort(GPP_B(:,2),'descend');
%     
% using the SVM regression
    [~,posSIF] = sort(SIF_Bsvm,'descend');
    [~,posGPP]  = sort(GPP_Bsvm,'descend');
        
        tickLabels={'None','Low','Moderate','High','Very High'};        
          for j=1:5
              tickLabels_ord{j} = tickLabels{posGPP(j)};
          end

          
          % Ordered using GPP
          f = figure;
          f.Position = [100 100 1500 600];
          for j=1:5
              yyaxis right
              h2 = bar(j,GPP_Bsvm(posGPP(j)),'FaceAlpha',0.5,'Facecolor',colors(posGPP(j),:));
              % Add the error bars for the GPP
              hold on
              errorbar(j,GPP_Bsvm(posGPP(j)),slopeGPP_SE{posGPP(j)},'k','LineStyle', 'none', 'LineWidth', 1.5);
              
              
              
              yyaxis left
              h1 = plot(j,SIF_Bsvm(posGPP(j)),'o','Markersize',20,'markerfacecolor',colors(posGPP(j),:),'markeredgecolor','k');
              hold on
              errorbar(j, SIF_Bsvm(posGPP(j)), slopeSIF_SE{posGPP(j)}, 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 10);
              
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
          xlim([0,6])
          title('GPP and cSIF anomalies response to air Temperature')
          lgnd = legend([h1,h2],'SIF','GPP','Location','east','orientation','horizontal');
          set(lgnd,'color',[0.8,0.8,0.8]);
          
          if save_flag
              switch years_
                  case 'initial'
                      saveas(gcf,[path_fig,'Thermokast_GPPcSIF_anom_changevsairT_order_2002_2016.png'])
                      saveas(gcf,[path_fig,'Thermokast_GPPcSIF_anom_changevsairT_order_2002_2016.fig'])
                      
                  case 'final'
                      saveas(gcf,[path_fig,'Thermokast_GPPcSIF_anom_changevsairT_order_2003_2018.png'])
                      saveas(gcf,[path_fig,'Thermokast_GPPcSIF_anom_changevsairT_order_2003_2018.fig'])
                      
                  case 'cSIF_all'
                      saveas(gcf,[path_fig,'Thermokast_GPPcSIF_anom_changevsairT_order_2003_2016.png'])
                      saveas(gcf,[path_fig,'Thermokast_GPPcSIF_anom_changevsairT_order_2003_2016.fig'])
                      
                  case 'cSIF_clear'
                      saveas(gcf,[path_fig,'Thermokast_GPPcSIFclear_anom_changevsairT_order_2003_2018.png'])
                      saveas(gcf,[path_fig,'Thermokast_GPPcSIFclear_anom_changevsairT_order_2003_2018.fig'])
              end
          end
% end

%% Plotting the slope values for SIF and GPP using different regressors 

f = figure;
f.Position = [100 100 1500 600];
for j=1:5
    h1 = bar((j-1)*4+1,SIF_B(j,2),'FaceAlpha',0.1,'Facecolor',colors(j,:));
    hold on
    h2 = bar((j-1)*4+2,SIF_Bsvm(j),'FaceAlpha',0.3,'Facecolor',colors(j,:));
    h3 = bar((j-1)*4+3,SIF_Blasso(j),'FaceAlpha',0.7,'Facecolor',colors(j,:));
    h4 = bar((j-1)*4+4,SIF_Benet(j),'FaceAlpha',1,'Facecolor',colors(j,:)); 
end
legend([h1,h2,h3,h4], 'OLS','SVM','LASSO','ENET')
xticks([2.5:4:5.*4])
xticklabels({'None','Low','Moderate','High','Very High'}); xtickangle(45)
set(get(gca, 'XAxis'), 'FontWeight', 'bold');
box on
grid on
grid minor
set(gca,'Fontsize',18)
% set(gca,'TickLength',[0.05, 0.005])
set(gca,'XMinorTick','on','YMinorTick','on')
title('SIF anomalies = \beta_0 + \beta_1 T anomalies')
ylabel('Slope linear regressors')


f = figure;
f.Position = [100 100 1500 600];
for j=1:5
    h1 = bar((j-1)*4+1,GPP_B(j,2),'FaceAlpha',0.1,'Facecolor',colors(j,:));
    hold on
    h2 = bar((j-1)*4+2,GPP_Bsvm(j),'FaceAlpha',0.3,'Facecolor',colors(j,:));
    h3 = bar((j-1)*4+3,GPP_Blasso(j),'FaceAlpha',0.7,'Facecolor',colors(j,:));
    h4 = bar((j-1)*4+4,GPP_Benet(j),'FaceAlpha',1,'Facecolor',colors(j,:)); 
end
legend([h1,h2,h3,h4], 'OLS','SVM','LASSO','ENET')
xticks([2.5:4:5.*4])
xticklabels({'None','Low','Moderate','High','Very High'}); xtickangle(45)
set(get(gca, 'XAxis'), 'FontWeight', 'bold');
box on
grid on
grid minor
set(gca,'Fontsize',18)
% set(gca,'TickLength',[0.05, 0.005])
set(gca,'XMinorTick','on','YMinorTick','on')
title('GPP anomalies = \beta_0 + \beta_1 T anomalies')
ylabel('Slope linear regressors')

if save_flag
    save([path_in,'Regressors_slope_thermokast.mat'],'SIF_B','SIF_Bsvm','SIF_Blasso','SIF_Benet',...
        'GPP_B','GPP_Bsvm','GPP_Blasso','GPP_Benet')
end


% --------------------------------------
% Paper figure 
% --------------------------------------


f = figure;
f.Position = [100 100 1000 500];
for j=1:5
       
    yyaxis right
    % GPP
    h1g = bar((j-1)*4+1,GPP_B(j,2),'FaceAlpha',0.1,'Facecolor',colors(j,:));
    hold on
    h2g =bar((j-1)*4+2,GPP_Bsvm(j),'FaceAlpha',0.3,'Facecolor',colors(j,:));
    h3g =bar((j-1)*4+3,GPP_Blasso(j),'FaceAlpha',0.7,'Facecolor',colors(j,:));
    h4g =bar((j-1)*4+4,GPP_Benet(j),'FaceAlpha',0.4,'Facecolor','w','Edgecolor',colors(j,:),'Linewidth',2);
%     hatchfill2(h4g,'cross','HatchAngle',30,'hatchcolor',colors(j,:),'FaceColor','none');
%     h4g.FaceColor = 'none';
    
    yyaxis left
    h1s = plot((j-1)*4+1,SIF_B(j,2),'o','Markersize',20,'markerfacecolor',colors(j,:),'markeredgecolor','k');
    hold on
    h2s =plot((j-1)*4+2,SIF_Bsvm(j),'o','Markersize',20,'markerfacecolor',colors(j,:),'markeredgecolor','k');
    h3s =plot((j-1)*4+3,SIF_Blasso(j),'o','Markersize',20,'markerfacecolor',colors(j,:),'markeredgecolor','k');
    h4s =plot((j-1)*4+4,SIF_Benet(j),'o','Markersize',20,'markerfacecolor',colors(j,:),'markeredgecolor','k');

end

xticks([2.5:4:5.*4])
xticklabels({'None','Low','Moderate','High','Very High'}); xtickangle(0)
set(get(gca, 'XAxis'), 'FontWeight', 'normal');
box on
grid on
grid minor
set(gca,'Fontsize',20)
% set(gca,'TickLength',[0.05, 0.005])
set(gca,'XMinorTick','on','YMinorTick','on')
title('SIF (GPP) anomalies = \beta_0 + \beta_1 T anomalies')
yyaxis right
ylabel('GPP anomaly vs. air T anomaly (\beta_1)')
yyaxis left
ylabel('SIF anomaly vs. air T anomaly (\beta_1)')
box on
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
yyaxis right 
ylim([0,0.3])
% Black bars to divide
yyaxis left
plot(4.5.*ones(size([3e-3:0.1e-3:11e-3])),[3e-3:0.1e-3:11e-3],'k-','Linewidth',4)
plot(8.5.*ones(size([3e-3:0.1e-3:11e-3])),[3e-3:0.1e-3:11e-3],'k-','Linewidth',4)
plot(12.5.*ones(size([3e-3:0.1e-3:11e-3])),[3e-3:0.1e-3:11e-3],'k-','Linewidth',4)
plot(16.5.*ones(size([3e-3:0.1e-3:11e-3])),[3e-3:0.1e-3:11e-3],'k-','Linewidth',4)
legend([h1s,h1g,h2g,h3g,h4g], 'SIF',...
    'OLS GPP','SVM GPP','LASSO GPP','ENET GPP')

if save_flag
    saveas(gcf,[path_fig,'Thermokast_GPPcSIFclear005_anom_changevsairT_ALLREG_2003_2018.png'])
    saveas(gcf,[path_fig,'Thermokast_GPPcSIFclear005_anom_changevsairT_ALLREG_2003_2018.fig'])
end

% --------------------------------------
% Paper figure combined
% --------------------------------------


f = figure;
f.Position = [100 100 1000 500];
for j=1:5
    
    % Compute the statistics
    mean_val_GPP = mean([GPP_B(j,2),GPP_Bsvm(j),GPP_Blasso(j),GPP_Benet(j)]);
    min_val_GPP  = min([GPP_B(j,2),GPP_Bsvm(j),GPP_Blasso(j),GPP_Benet(j)]);
    max_val_GPP  = max([GPP_B(j,2),GPP_Bsvm(j),GPP_Blasso(j),GPP_Benet(j)]);
    diff_min_GPP = abs(mean_val_GPP-min_val_GPP);
    diff_max_GPP = abs(max_val_GPP-mean_val_GPP);
        
    
    mean_val_SIF = mean([SIF_B(j,2),SIF_Bsvm(j),SIF_Blasso(j),SIF_Benet(j)]);    
    min_val_SIF  = min([SIF_B(j,2),SIF_Bsvm(j),SIF_Blasso(j),SIF_Benet(j)]);
    max_val_SIF  = max([SIF_B(j,2),SIF_Bsvm(j),SIF_Blasso(j),SIF_Benet(j)]);
    diff_min_SIF = abs(mean_val_SIF-min_val_SIF);
    diff_max_SIF = abs(max_val_SIF-mean_val_SIF);
    
    
    
    yyaxis right
    % GPP
    h1g = plot(j,mean_val_GPP,'s','Markersize',20,'markerfacecolor','k','markeredgecolor','k');
%     h1g = bar(j,mean_val_GPP,'FaceAlpha',0.5,'Facecolor',colors(j,:));
    hold on
    errorbar(j,mean_val_GPP, diff_min_GPP, diff_max_GPP, 'r', 'linestyle', 'none');
    hold on
    ylim([0.1,0.3])
    
    yyaxis left
    h1s = plot(j,mean_val_SIF,'o','Markersize',20,'markerfacecolor',colors(j,:),'markeredgecolor','k');
    hold on 
    errorbar(j, mean_val_SIF, diff_min_SIF, diff_max_SIF, 'k--', 'linestyle', 'none');
    ylim([0,11e-3])
    
end
xlim([0,6])
xticks([1:5])
xticklabels({'None','Low','Moderate','High','Very High'}); xtickangle(0)
set(get(gca, 'XAxis'), 'FontWeight', 'normal');
box on
grid on
grid minor
set(gca,'Fontsize',20)
% set(gca,'TickLength',[0.05, 0.005])
set(gca,'XMinorTick','on','YMinorTick','on')
title('SIF (GPP) anomalies = \beta_0 + \beta_1 T anomalies')
yyaxis right
ylabel('GPP anomaly vs. air T anomaly (\beta_1)')
yyaxis left
ylabel('SIF anomaly vs. air T anomaly (\beta_1)')
box on
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';


if save_flag
    saveas(gcf,[path_fig,'Thermokast_GPPcSIFclear005_anom_changevsairT_ALLREG_2003_2018.png'])
    saveas(gcf,[path_fig,'Thermokast_GPPcSIFclear005_anom_changevsairT_ALLREG_2003_2018.fig'])
end


% ------------------------------
% Alternative
%% ------------------------------


f = figure;
f.Position = [100 100 700 450];
for j=1:5
    
    % Compute the statistics for GPP
    mean_val_GPP = mean([GPP_B(j,2),GPP_Bsvm(j),GPP_Blasso(j),GPP_Benet(j)]);
    min_val_GPP  = min([GPP_B(j,2),GPP_Bsvm(j),GPP_Blasso(j),GPP_Benet(j)]);
    max_val_GPP  = max([GPP_B(j,2),GPP_Bsvm(j),GPP_Blasso(j),GPP_Benet(j)]);
    diff_min_GPP = abs(mean_val_GPP-min_val_GPP);
    diff_max_GPP = abs(max_val_GPP-mean_val_GPP);
        
    
    % Compute the statistics for SIF
    mean_val_SIF = mean([SIF_B(j,2),SIF_Bsvm(j),SIF_Blasso(j),SIF_Benet(j)]);    
    min_val_SIF  = min([SIF_B(j,2),SIF_Bsvm(j),SIF_Blasso(j),SIF_Benet(j)]);
    max_val_SIF  = max([SIF_B(j,2),SIF_Bsvm(j),SIF_Blasso(j),SIF_Benet(j)]);
    diff_min_SIF = abs(mean_val_SIF-min_val_SIF);
    diff_max_SIF = abs(max_val_SIF-mean_val_SIF);
    
    
    % Create shaded region for GPP (right y-axis)
    yyaxis right
    hold on;
    if j == 2
        fill([j-0.1 j+0.1 j+0.1 j-0.1], [min_val_GPP min_val_GPP max_val_GPP max_val_GPP], colors(j,:), 'FaceAlpha', 0.4, 'EdgeColor', 'k','Linestyle',':','Marker','none');
        hg = plot(j, mean_val_GPP, 's', 'Markersize', 20, 'markerfacecolor', colors(j,:), 'markeredgecolor', 'k');
    else
        fill([j-0.1 j+0.1 j+0.1 j-0.1], [min_val_GPP min_val_GPP max_val_GPP max_val_GPP], colors(j,:), 'FaceAlpha', 0.4, 'EdgeColor', 'k','Linestyle',':','Marker','none');
        plot(j, mean_val_GPP, 's', 'Markersize', 20, 'markerfacecolor', colors(j,:), 'markeredgecolor', 'k');
    end
    ylim([0.1,0.3])
    
    % Create shaded region for SIF (left y-axis)
    yyaxis left
    hold on;
    if j==2
        hh = fill([j-0.1 j+0.1 j+0.1 j-0.1], [min_val_SIF min_val_SIF max_val_SIF max_val_SIF], colors(j,:), 'FaceAlpha', 0.4, 'EdgeColor', 'k','Linestyle',':','Marker','none');
        hs = plot(j, mean_val_SIF, 'o', 'Markersize', 20, 'markerfacecolor', colors(j,:), 'markeredgecolor', 'k');
    else
        fill([j-0.1 j+0.1 j+0.1 j-0.1], [min_val_SIF min_val_SIF max_val_SIF max_val_SIF], colors(j,:), 'FaceAlpha', 0.4, 'EdgeColor', 'k','Linestyle',':','Marker','none');
        plot(j, mean_val_SIF, 'o', 'Markersize', 20, 'markerfacecolor', colors(j,:), 'markeredgecolor', 'k');
    end
    ylim([0,11e-3])
    
end

% Set x-axis limits and labels
xlim([0.5,5.5])
xticks(1:5)
xticklabels({'None','Low','Moderate','High','Very High'});
xtickangle(0)
set(get(gca, 'XAxis'), 'FontWeight', 'normal');

% Formatting
box on
set(gca,'Fontsize',20)

title('SIF (GPP) anomalies = \beta_0 + \beta_1 T anomalies')

% Y-axis labels
yyaxis right
ylabel('GPP anomaly vs. air T anomaly (\beta_1)')
set(gca,'XMinorTick','off','YMinorTick','on')
yyaxis left
ylabel('SIF anomaly vs. air T anomaly (\beta_1)')
set(gca,'XMinorTick','off','YMinorTick','on')

% Set axis colors
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';



yyaxis left
% plot(1.5.*ones(size([0:0.1e-3:11e-3])),[0:0.1e-3:11e-3],'k-','Linewidth',2)
% plot(2.5.*ones(size([0:0.1e-3:11e-3])),[0:0.1e-3:11e-3],'k-','Linewidth',2)
% plot(3.5.*ones(size([0:0.1e-3:11e-3])),[0:0.1e-3:11e-3],'k-','Linewidth',2)
% plot(4.5.*ones(size([0:0.1e-3:11e-3])),[0:0.1e-3:11e-3],'k-','Linewidth',2)

% legend 
legend([hs,hg,hh], 'SIF \beta_1','GPP \beta_1','min-max','location','Northwest')
legend boxoff;
set(gca,'XMinorTick','off','YMinorTick','on')


if save_flag
    saveas(gcf,[path_fig,'Thermokast_GPPcSIFclear005_anom_changevsairT_ALLREG_Combined_2003_2018.png'])
    saveas(gcf,[path_fig,'Thermokast_GPPcSIFclear005_anom_changevsairT_ALLREG_Combined_2003_2018.fig'])
end






















%% Plotting ordering by method but ALL the regressors

method_name = 'OLS';  %'OLS','SVM','LASSO','BENET'

% if plotting == 1
switch method_name
    case 'OLS'
        % using the OLS regression
        [~,posSIF] = sort(SIF_B(:,2),'descend');
        [~,posGPP]  = sort(GPP_B(:,2),'descend');
    case 'SVM'        
        % using the SVM regression
        [~,posSIF] = sort(SIF_Bsvm,'descend');
        [~,posGPP]  = sort(GPP_Bsvm,'descend');
    case 'LASSO'   
        [~,posSIF] = sort(SIF_Blasso,'descend');
        [~,posGPP]  = sort(GPP_Blasso,'descend');
    case 'BENET'
        [~,posSIF] = sort(SIF_Benet,'descend');
        [~,posGPP]  = sort(GPP_Benet,'descend');
end


f = figure;
f.Position = [100 100 1500 600];
for j=1:5
       
    yyaxis right
    % GPP
    h1g = bar((j-1)*4+1,GPP_B(posGPP(j),2),'FaceAlpha',0.1,'Facecolor',colors(posGPP(j),:));
    hold on
    h2g =bar((j-1)*4+2,GPP_Bsvm(posGPP(j)),'FaceAlpha',0.3,'Facecolor',colors(posGPP(j),:));
    h3g =bar((j-1)*4+3,GPP_Blasso(posGPP(j)),'FaceAlpha',0.7,'Facecolor',colors(posGPP(j),:));
    h4g =bar((j-1)*4+4,GPP_Benet(posGPP(j)),'FaceAlpha',0.4,'Facecolor','w','Edgecolor',colors(posGPP(j),:),'Linewidth',2);
%     hatchfill2(h4g,'cross','HatchAngle',30,'hatchcolor',colors(j,:),'FaceColor','none');
%     h4g.FaceColor = 'none';
    
    yyaxis left
    h1s = plot((j-1)*4+1,SIF_B(posGPP(j),2),'o','Markersize',20,'markerfacecolor',colors(posGPP(j),:),'markeredgecolor','k');
    hold on
    h2s =plot((j-1)*4+2,SIF_Bsvm(posGPP(j)),'o','Markersize',20,'markerfacecolor',colors(posGPP(j),:),'markeredgecolor','k');
    h3s =plot((j-1)*4+3,SIF_Blasso(posGPP(j)),'o','Markersize',20,'markerfacecolor',colors(posGPP(j),:),'markeredgecolor','k');
    h4s =plot((j-1)*4+4,SIF_Benet(posGPP(j)),'o','Markersize',20,'markerfacecolor',colors(posGPP(j),:),'markeredgecolor','k');

end
legend([h1s,h1g,h2g,h3g,h4g], 'SIF',...
    'OLS GPP','SVM GPP','LASSO GPP','BENET GPP')
xticks([2.5:4:5.*4])
for j=1:5
    tickLabels_ord{j} = tickLabels{posGPP(j)};
end
xticklabels(tickLabels_ord)
xtickangle(45)
set(get(gca, 'XAxis'), 'FontWeight', 'bold');
box on
grid on
grid minor
set(gca,'Fontsize',18)
% set(gca,'TickLength',[0.05, 0.005])
set(gca,'XMinorTick','on','YMinorTick','on')
title('SIF (GPP) anomalies = \beta_0 + \beta_1 T anomalies')
yyaxis right
ylabel('GPP anomaly vs. air T anomaly (\beta_1)')
yyaxis left
ylabel('SIF anomaly vs. air T anomaly (\beta_1)')
box on
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
yyaxis right 
ylim([0,0.3])

yyaxis left
plot(4.5.*ones(size([3e-3:0.1e-3:11e-3])),[3e-3:0.1e-3:11e-3],'k-','Linewidth',4)
plot(8.5.*ones(size([3e-3:0.1e-3:11e-3])),[3e-3:0.1e-3:11e-3],'k-','Linewidth',4)
plot(12.5.*ones(size([3e-3:0.1e-3:11e-3])),[3e-3:0.1e-3:11e-3],'k-','Linewidth',4)
plot(16.5.*ones(size([3e-3:0.1e-3:11e-3])),[3e-3:0.1e-3:11e-3],'k-','Linewidth',4)
legend([h1s,h1g,h2g,h3g,h4g], 'SIF',...
    'OLS GPP','SVM GPP','LASSO GPP','BENET GPP')

if save_flag
    saveas(gcf,[path_fig,'Thermokast_GPPcSIFclear005_ord_anom_changevsairT_ALLREG_2003_2018.png'])
    saveas(gcf,[path_fig,'Thermokast_GPPcSIFclear005_ord_anom_changevsairT_ALLREG_2003_2018.fig'])
end

%% Plotting ordering by method

method_name = 'LASSO';  %'OLS','SVM','LASSO','BENET'

% if plotting == 1
switch method_name
    case 'OLS'
        % using the OLS regression
        [~,posSIF] = sort(SIF_B(:,2),'descend');
        [~,posGPP]  = sort(GPP_B(:,2),'descend');
    case 'SVM'        
        % using the SVM regression
        [~,posSIF] = sort(SIF_Bsvm,'descend');
        [~,posGPP]  = sort(GPP_Bsvm,'descend');
    case 'LASSO'   
        [~,posSIF] = sort(SIF_Blasso,'descend');
        [~,posGPP]  = sort(GPP_Blasso,'descend');
    case 'BENET'
        [~,posSIF] = sort(SIF_Benet,'descend');
        [~,posGPP]  = sort(GPP_Benet,'descend');
end

    tickLabels={'None','Low','Moderate','High','Very High'};        
        
          for j=1:5
              tickLabels_ord{j} = tickLabels{posGPP(j)};
          end

          
          % Ordered using GPP
          f = figure;
          f.Position = [100 100 1500 600];
          for j=1:5
              
              switch method_name
                  case 'OLS'
                      yyaxis right
                      h2 = bar(j,GPP_B(posGPP(j),2),'FaceAlpha',0.5,'Facecolor',colors(posGPP(j),:));
                      hold on
                      yyaxis left
                      h1 = plot(j,SIF_B(posGPP(j),2),'o','Markersize',20,'markerfacecolor',colors(posGPP(j),:),'markeredgecolor','k');
                      hold on
                  case 'SVM'
                      yyaxis right
                      h2 = bar(j,GPP_Bsvm(posGPP(j)),'FaceAlpha',0.5,'Facecolor',colors(posGPP(j),:));
                      hold on
                      yyaxis left
                      h1 = plot(j,SIF_Bsvm(posGPP(j)),'o','Markersize',20,'markerfacecolor',colors(posGPP(j),:),'markeredgecolor','k');
                      hold on
                  case 'LASSO'
                      yyaxis right
                      h2 = bar(j,GPP_Blasso(posGPP(j)),'FaceAlpha',0.5,'Facecolor',colors(posGPP(j),:));
                      hold on
                      yyaxis left
                      h1 = plot(j,SIF_Blasso(posGPP(j)),'o','Markersize',20,'markerfacecolor',colors(posGPP(j),:),'markeredgecolor','k');
                      hold on  
                  case 'BENET'
                      yyaxis right
                      h2 = bar(j,GPP_Benet(posGPP(j)),'FaceAlpha',0.5,'Facecolor',colors(posGPP(j),:));
                      hold on
                      yyaxis left
                      h1 = plot(j,SIF_Benet(posGPP(j)),'o','Markersize',20,'markerfacecolor',colors(posGPP(j),:),'markeredgecolor','k');
                      hold on  
              end
              
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
          xticks([1:5])
          xticklabels(tickLabels_ord)
          xtickangle(45)
          grid minor
          set(gca,'Fontsize',18)
          set(gca,'TickLength',[0.05, 0.005])
          set(gca,'XMinorTick','on','YMinorTick','on')
          xlim([0,16])
          title(['GPP and cSIF anomalies response to air Temperature  [',method_name,']'])
          lgnd = legend([h1,h2],'SIF','GPP','Location','east','orientation','horizontal');
          set(lgnd,'color',[0.8,0.8,0.8]);
          
          if save_flag
              saveas(gcf,[path_fig,'Thermokast_GPPcSIFclear_anom_changevsairT_order_2003_2018',method_name,'.png'])
              saveas(gcf,[path_fig,'Thermokast_GPPcSIFclear_anom_changevsairT_order_2003_2018',method_name,'.fig'])
          end