% ###################################################################
% ####### plotPolarPlotAnomCSIFAll_class   ####
% This script reads the cSIF database anomalies from 2003-2016 (or 2003-2018) and computes 
% the OLS regression with the corresponding airT anomalies for the same
% years. 
% Data was extracted from the cSIF database generated with OCO-2 SIF data and 
% MODIS surface reflectance every 4-day at 0.05 degrees. 
% This database also disentangles the sunlit and all sky SIF prediction
% 
% %Zhang, Y., Joiner, J., Alemohammad, S. H., Zhou, S., & Gentine, P. (2018).
% A global spatially contiguous solar-induced fluorescence (CSIF) dataset 
% using neural networks. Biogeosciences, 15(19), 5779-5800.
% https://bg.copernicus.org/articles/15/5779/2018/

% Author: Neus Sabater
% e-mail: neus.sabater@fmi.fi 
% version 0: (19/02/2024) 
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

% classification in BAWLD
CLASS_int = ncread([path_in,'BAWLD_dataset.nc'],'CLASS_int');


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
    f = figure;
    f.Position = [100 100 2000 1200];
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
        
        % matlab OLS linear regressor
        b = X\y;
        GPP_B(j,:) = b;
        yCalc2 = X*b;
        Rsq2    = 1 - sum((y - yCalc2).^2)/sum((y - mean(y)).^2);
        
        % Compute other regressors info
        mdl   = fitlm(x,y); % OLS linear
        disp('Report statistics GPP regressor')
        % Store the model and its statistics in the structure
        resultsGPP(j).Model           = mdl;  % Store the entire model object
        resultsGPP(j).Coefficients    = mdl.Coefficients;  % Coefficients and related statistics
        resultsGPP(j).R2              = mdl.Rsquared.Ordinary;  % R^2
        resultsGPP(j).AdjustedR2      = mdl.Rsquared.Adjusted;  % Adjusted R^2
        resultsGPP(j).RMSE            = mdl.RMSE;  % Root Mean Squared Error
        resultsGPP(j).PValues         = mdl.Coefficients.pValue;  % p-values for coefficients
        resultsGPP(j).NumObservations = mdl.NumObservations;  % Number of observations
        resultsGPP(j).ErrorDOF        = mdl.DFE;  % Error degrees of freedom
        resultsGPP(j).beta1           = mdl.Coefficients.Estimate(2); % Error associated with the slope
        aux                           = mdl.coefCI;
        resultsGPP(j).Conf_int        = aux(2,:); % Confidence interval associated with the slope
        resultsGPP(j).SE_int          = mdl.Coefficients.SE(2); % Standard error associated with the slope
        resultsGPP(j).pValue_x1       = mdl.Coefficients.pValue(2);
                
    
        
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
                saveas(gcf,[path_fig,'BAWLD_GPP_class_scatter_2002_2016.png'])
                saveas(gcf,[path_fig,'BAWLD_GPP_class_scatter_2002_2016.fig'])
                
            case 'final'
                saveas(gcf,[path_fig,'BAWLD_GPP_class_scatter_2003_2018.png'])
                saveas(gcf,[path_fig,'BAWLD_GPP_class_scatter_2003_2018.fig'])
                
            case 'cSIF_all'
                saveas(gcf,[path_fig,'BAWLD_GPP_class_scatter_2003_2016.png'])
                saveas(gcf,[path_fig,'BAWLD_GPP_class_scatter_2003_2016.fig'])
                
            case 'cSIF_clear'
                saveas(gcf,[path_fig,'BAWLD_GPPclear_class_scatter_2003_2018.png'])
                saveas(gcf,[path_fig,'BAWLD_GPPclear_class_scatter_2003_2018.fig'])
                
        end
    end
    
    
%   cSIF anomalies 
    f = figure;
    f.Position = [100 100 2000 1200];
    for j=1:15
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
        
        % Compute other regressors info
        mdl   = fitlm(x,y); % OLS linear
        disp('Report statistics cSIF')
        % Store the model and its statistics in the structure
        resultsSIF(j).Model           = mdl;  % Store the entire model object
        resultsSIF(j).Coefficients    = mdl.Coefficients;  % Coefficients and related statistics
        resultsSIF(j).R2              = mdl.Rsquared.Ordinary;  % R^2
        resultsSIF(j).AdjustedR2      = mdl.Rsquared.Adjusted;  % Adjusted R^2
        resultsSIF(j).RMSE            = mdl.RMSE;  % Root Mean Squared Error
        resultsSIF(j).PValues         = mdl.Coefficients.pValue;  % p-values for coefficients
        resultsSIF(j).NumObservations = mdl.NumObservations;  % Number of observations
        resultsSIF(j).ErrorDOF        = mdl.DFE;  % Error degrees of freedom
        resultsSIF(j).beta1           = mdl.Coefficients.Estimate(2); % Error associated with the slope
        aux                           = mdl.coefCI;
        resultsSIF(j).Conf_int        = aux(2,:); % Confidence interval associated with the slope
        resultsSIF(j).SE_int          = mdl.Coefficients.SE(2); % Standard error associated with the slope
        resultsSIF(j).pValue_x1       = mdl.Coefficients.pValue(2);

       
        
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
                saveas(gcf,[path_fig,'BAWLD_cSIF_class_scatter_2002_2016.png'])
                saveas(gcf,[path_fig,'BAWLD_cSIF_class_scatter_2002_2016.fig'])
                
            case 'final'
                saveas(gcf,[path_fig,'BAWLD_cSIF_class_scatter_2003_2018.png'])
                saveas(gcf,[path_fig,'BAWLD_cSIF_class_scatter_2003_2018.fig'])
                
            case 'cSIF_all'
                saveas(gcf,[path_fig,'BAWLD_cSIF_class_scatter_2003_2016.png'])
                saveas(gcf,[path_fig,'BAWLD_cSIF_class_scatter_2003_2016.fig'])
                
            case 'cSIF_clear'
                saveas(gcf,[path_fig,'BAWLD_cSIFclear_class_scatter_2003_2018.png'])
                saveas(gcf,[path_fig,'BAWLD_cSIFclear_class_scatter_2003_2018.fig'])
        end
    end
    
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
                saveas(gcf,[path_fig,'BAWLD_cSIF_linearreg_2003_2018.png'])
                saveas(gcf,[path_fig,'BAWLD_cSIF_linearreg_2003_2018.fig'])
                
            case 'cSIF_all'
                saveas(gcf,[path_fig,'BAWLD_cSIF_linearreg_2003_2016.png'])
                saveas(gcf,[path_fig,'BAWLD_cSIF_linearreg_2003_2016.fig'])
                
            case 'cSIF_clear'
                saveas(gcf,[path_fig,'BAWLD_cSIFclear_linearreg_2003_2018.png'])
                saveas(gcf,[path_fig,'BAWLD_cSIFclear_linearreg_2003_2018.fig'])
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
                
            case 'cSIF_all'
                saveas(gcf,[path_fig,'BAWLD_GPP_linearreg_2003_2016.png'])
                saveas(gcf,[path_fig,'BAWLD_GPP_linearreg_2003_2016.fig'])
                
            case 'cSIF_clear'
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
                saveas(gcf,[path_fig,'BAWLD_GPPcSIF_anom_changevsairT_2002_2016.png'])
                saveas(gcf,[path_fig,'BAWLD_GPPcSIF_anom_changevsairT_2002_2016.fig'])
                
            case 'final'
                saveas(gcf,[path_fig,'BAWLD_GPPcSIF_anom_changevsairT_2003_2018.png'])
                saveas(gcf,[path_fig,'BAWLD_GPPcSIF_anom_changevsairT_2003_2018.fig'])
            
            case 'cSIF_all'
                saveas(gcf,[path_fig,'BAWLD_GPPcSIF_anom_changevsairT_2003_2016.png'])
                saveas(gcf,[path_fig,'BAWLD_GPPcSIF_anom_changevsairT_2003_2016.fig'])
                
            case 'cSIF_clear'
                saveas(gcf,[path_fig,'BAWLD_GPPcSIFclear_anom_changevsairT_2003_2018.png'])
                saveas(gcf,[path_fig,'BAWLD_GPPcSIFclear_anom_changevsairT_2003_2018.fig'])
                
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
%             h2 = bar(j,GPP_B(posGPP(j),2),'FaceAlpha',0.5,'Facecolor',colors(posGPP(j),:));
              h2 = bar(j,GPP_Bsvm(posGPP(j)),'FaceAlpha',0.5,'Facecolor',colors(posGPP(j),:));

              hold on
              yyaxis left
%               h1 = plot(j,SIF_B(posGPP(j),2),'o','Markersize',20,'markerfacecolor',colors(posGPP(j),:),'markeredgecolor','k');
              h1 = plot(j,SIF_Bsvm(posGPP(j)),'o','Markersize',20,'markerfacecolor',colors(posGPP(j),:),'markeredgecolor','k');

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
          title('GPP and cSIF anomalies response to air Temperature')
          lgnd = legend([h1,h2],'SIF','GPP','Location','east','orientation','horizontal');
          set(lgnd,'color',[0.8,0.8,0.8]);
          
          if save_flag
              switch years_
                  case 'initial'
                      saveas(gcf,[path_fig,'BAWLD_GPPcSIF_anom_changevsairT_order_2002_2016.png'])
                      saveas(gcf,[path_fig,'BAWLD_GPPcSIF_anom_changevsairT_order_2002_2016.fig'])
                      
                  case 'final'
                      saveas(gcf,[path_fig,'BAWLD_GPPcSIF_anom_changevsairT_order_2003_2018.png'])
                      saveas(gcf,[path_fig,'BAWLD_GPPcSIF_anom_changevsairT_order_2003_2018.fig'])
                      
                  case 'cSIF_all'
                      saveas(gcf,[path_fig,'BAWLD_GPPcSIF_anom_changevsairT_order_2003_2016.png'])
                      saveas(gcf,[path_fig,'BAWLD_GPPcSIF_anom_changevsairT_order_2003_2016.fig'])
                      
                  case 'cSIF_clear'
                      saveas(gcf,[path_fig,'BAWLD_GPPcSIFclear_anom_changevsairT_order_2003_2018.png'])
                      saveas(gcf,[path_fig,'BAWLD_GPPcSIFclear_anom_changevsairT_order_2003_2018.fig'])
              end
          end
% end

%% Plotting the slope values for SIF and GPP using different regressors 

f = figure;
f.Position = [100 100 1500 600];
for j=1:15
    h1 = bar((j-1)*4+1,SIF_B(j,2),'FaceAlpha',0.1,'Facecolor',colors(j,:));
    hold on
    h2 = bar((j-1)*4+2,SIF_Bsvm(j),'FaceAlpha',0.3,'Facecolor',colors(j,:));
    h3 = bar((j-1)*4+3,SIF_Blasso(j),'FaceAlpha',0.7,'Facecolor',colors(j,:));
    h4 = bar((j-1)*4+4,SIF_Benet(j),'FaceAlpha',1,'Facecolor',colors(j,:)); 
end
legend([h1,h2,h3,h4], 'OLS','SVM','LASSO','ENET')
xticks([2.5:4:15.*4])
xticklabels({'Permafrost Peatlands','Sparse Boreal Peatlands','Rivers',...
    'Glaciers','Upland Tundra','Common Boreal Peatlands','Large Lakes',...
    'Lake-rich wetlands','Dominant Boreal Peatlands','Wetland-rich Tundra',...
    'Alpine and Tundra Barrens','Wetland and Lake-rich Tundra','Lake-rich Shield','Upland Boreal','Wetland and Lake-rich Yedoma Tundra'});
xtickangle(45)
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
for j=1:15
    h1 = bar((j-1)*4+1,GPP_B(j,2),'FaceAlpha',0.1,'Facecolor',colors(j,:));
    hold on
    h2 = bar((j-1)*4+2,GPP_Bsvm(j),'FaceAlpha',0.3,'Facecolor',colors(j,:));
    h3 = bar((j-1)*4+3,GPP_Blasso(j),'FaceAlpha',0.7,'Facecolor',colors(j,:));
    h4 = bar((j-1)*4+4,GPP_Benet(j),'FaceAlpha',1,'Facecolor',colors(j,:)); 
end
legend([h1,h2,h3,h4], 'OLS','SVM','LASSO','ENET')
xticks([2.5:4:15.*4])
xticklabels({'Permafrost Peatlands','Sparse Boreal Peatlands','Rivers',...
    'Glaciers','Upland Tundra','Common Boreal Peatlands','Large Lakes',...
    'Lake-rich wetlands','Dominant Boreal Peatlands','Wetland-rich Tundra',...
    'Alpine and Tundra Barrens','Wetland and Lake-rich Tundra','Lake-rich Shield','Upland Boreal','Wetland and Lake-rich Yedoma Tundra'});
xtickangle(45)
set(get(gca, 'XAxis'), 'FontWeight', 'bold');
box on
grid on
grid minor
set(gca,'Fontsize',18)
% set(gca,'TickLength',[0.05, 0.005])
set(gca,'XMinorTick','on','YMinorTick','on')
title('GPP anomalies = \beta_0 + \beta_1 T anomalies')
ylabel('Slope linear regressors')


save([path_in,'Regressors_slope.mat'],'SIF_B','SIF_Bsvm','SIF_Blasso','SIF_Benet',...
    'GPP_B','GPP_Bsvm','GPP_Blasso','GPP_Benet')



f = figure;
f.Position = [100 100 1500 600];
for j=1:15
       
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
legend([h1s,h1g,h2g,h3g,h4g], 'SIF',...
    'OLS GPP','SVM GPP','LASSO GPP','ENET GPP')
xticks([2.5:4:15.*4])
xticklabels({'Permafrost Peatlands','Sparse Boreal Peatlands','Rivers',...
    'Glaciers','Upland Tundra','Common Boreal Peatlands','Large Lakes',...
    'Lake-rich wetlands','Dominant Boreal Peatlands','Wetland-rich Tundra',...
    'Alpine and Tundra Barrens','Wetland and Lake-rich Tundra','Lake-rich Shield','Upland Boreal','Wetland and Lake-rich Yedoma Tundra'});
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

saveas(gcf,[path_fig,'BAWLD_GPPcSIFclear005_anom_changevsairT_ALLREG_2003_2018.png'])
saveas(gcf,[path_fig,'BAWLD_GPPcSIFclear005_anom_changevsairT_ALLREG_2003_2018.fig'])
%% Plotting ordering by method but ALL the regressors

method_name = 'OLS';  %'OLS','SVM','LASSO','ENET'

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
for j=1:15
       
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
    'OLS GPP','SVM GPP','LASSO GPP','ENET GPP')
xticks([2.5:4:15.*4])
for j=1:15
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

saveas(gcf,[path_fig,'BAWLD_GPPcSIFclear005_ord_anom_changevsairT_ALLREG_2003_2018.png'])
saveas(gcf,[path_fig,'BAWLD_GPPcSIFclear005_ord_anom_changevsairT_ALLREG_2003_2018.fig'])


%% Plotting ordering by method

method_name = 'LASSO';  %'OLS','SVM','LASSO','ENET'

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
          xticks([1:15])
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
          
          
          saveas(gcf,[path_fig,'BAWLD_GPPcSIFclear_anom_changevsairT_order_2003_2018',method_name,'.png'])
          saveas(gcf,[path_fig,'BAWLD_GPPcSIFclear_anom_changevsairT_order_2003_2018',method_name,'.fig'])


%% Paper figure combined 
f = figure;
f.Position = [100 100 900 450];

method_name = 'OLS';  %'OLS','SVM','LASSO','ENET'

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

for j=1:15
    
    % Compute the statistics for GPP
    mean_val_GPP = mean([GPP_B(posGPP(j),2),GPP_Bsvm(posGPP(j)),GPP_Blasso(posGPP(j)),GPP_Benet(posGPP(j))]);
    min_val_GPP  = min([GPP_B(posGPP(j),2),GPP_Bsvm(posGPP(j)),GPP_Blasso(posGPP(j)),GPP_Benet(posGPP(j))]);
    max_val_GPP  = max([GPP_B(posGPP(j),2),GPP_Bsvm(posGPP(j)),GPP_Blasso(posGPP(j)),GPP_Benet(posGPP(j))]);
    diff_min_GPP = abs(mean_val_GPP-min_val_GPP);
    diff_max_GPP = abs(max_val_GPP-mean_val_GPP);
        
    
    % Compute the statistics for SIF
    mean_val_SIF = mean([SIF_B(posGPP(j),2),SIF_Bsvm(posGPP(j)),SIF_Blasso(posGPP(j)),SIF_Benet(posGPP(j))]);    
    min_val_SIF  = min([SIF_B(posGPP(j),2),SIF_Bsvm(posGPP(j)),SIF_Blasso(posGPP(j)),SIF_Benet(posGPP(j))]);
    max_val_SIF  = max([SIF_B(posGPP(j),2),SIF_Bsvm(posGPP(j)),SIF_Blasso(posGPP(j)),SIF_Benet(posGPP(j))]);
    diff_min_SIF = abs(mean_val_SIF-min_val_SIF);
    diff_max_SIF = abs(max_val_SIF-mean_val_SIF);
    
    
    % Create shaded region for GPP (right y-axis)
    yyaxis right
    hold on;
    if j == 2
        fill([j-0.1 j+0.1 j+0.1 j-0.1], [min_val_GPP min_val_GPP max_val_GPP max_val_GPP], colors(posGPP(j),:), 'FaceAlpha', 0.4, 'EdgeColor', 'k','Linestyle',':','Marker','none');
        hg = plot(j, mean_val_GPP, 's', 'Markersize', 20, 'markerfacecolor', colors(posGPP(j),:), 'markeredgecolor', 'k');
    else
        fill([j-0.1 j+0.1 j+0.1 j-0.1], [min_val_GPP min_val_GPP max_val_GPP max_val_GPP], colors(posGPP(j),:), 'FaceAlpha', 0.4, 'EdgeColor', 'k','Linestyle',':','Marker','none');
        plot(j, mean_val_GPP, 's', 'Markersize', 20, 'markerfacecolor', colors(posGPP(j),:), 'markeredgecolor', 'k');
    end
    ylim([0.1,0.3])
    
    % Create shaded region for SIF (left y-axis)
    yyaxis left
    hold on;
    if j==2
        hh = fill([j-0.1 j+0.1 j+0.1 j-0.1], [min_val_SIF min_val_SIF max_val_SIF max_val_SIF], colors(posGPP(j),:), 'FaceAlpha', 0.4, 'EdgeColor', 'k','Linestyle',':','Marker','none');
        hs = plot(j, mean_val_SIF, 'o', 'Markersize', 20, 'markerfacecolor', colors(posGPP(j),:), 'markeredgecolor', 'k');
    else
        fill([j-0.1 j+0.1 j+0.1 j-0.1], [min_val_SIF min_val_SIF max_val_SIF max_val_SIF], colors(posGPP(j),:), 'FaceAlpha', 0.4, 'EdgeColor', 'k','Linestyle',':','Marker','none');
        plot(j, mean_val_SIF, 'o', 'Markersize', 20, 'markerfacecolor', colors(posGPP(j),:), 'markeredgecolor', 'k');
    end
    
end

% Set x-axis limits and labels
yyaxis right 
ylim([0,0.3])
yyaxis left
ylim([0,25e-3])
xticks(1:15)
tickLabels = {'Permafrost Peatlands','Sparse Boreal Peatlands','Rivers',...
    'Glaciers','Upland Tundra','Common Boreal Peatlands','Large Lakes',...
    'Lake-rich wetlands','Dominant Boreal Peatlands','Wetland-rich Tundra',...
    'Alpine and Tundra Barrens','Wetland and Lake-rich Tundra','Lake-rich Shield','Upland Boreal','Wetland and Lake-rich Yedoma Tundra'};

for j=1:15
    tickLabels_ord{j} = tickLabels{posGPP(j)};
end
xticklabels(tickLabels_ord)
set(get(gca, 'XAxis'), 'FontWeight', 'bold');
xtickangle(45)
grid minor
set(gca,'Fontsize',18)
set(gca,'TickLength',[0.01, 0.005])
set(gca,'XMinorTick','on','YMinorTick','on')
xlim([0,16])
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

lgnd = legend([hg,hs,hh],'cSIF (clear) \beta_1','GPP \beta_1','min-max','Location','NorthEast','orientation','vertical');
set(lgnd,'color',[1,1,1],'EdgeColor',[0,0,0]);


saveas(gcf,[path_fig,'BAWLD_GPPcSIFclear_anom_changevsairT_Combined_order_2003_2018',method_name,'.png'])
saveas(gcf,[path_fig,'BAWLD_GPPcSIFclear_anom_changevsairT_Combined_order_2003_2018',method_name,'.fig'])
