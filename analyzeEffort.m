% Analyze sourcedata (.xlsx) from Effort Task. Goal is to quantify how
% willingness to exert effort in social and monetary conditions is
% influenced by value.

% Notes:
% domain is the same thing as "condition"
% 1020, 1021, 1023 <- social subjects on different reward scale.
% need one subject list


clear; close all;
maindir = pwd;
warning off all

% loop through all participants
subs = {''};

% first loop through monetary domain, then social
domains = {'monetary', 'social'};
for d = 1:length(domains)
    domain = domains{d};
    
    % build path for data and create list of files
    sourcedatadir = fullfile(maindir,'data',domain);
    sourcedata = dir([sourcedatadir '/*.xls*']);
    sourcedata = struct2cell(sourcedata);
    sourcedata = sourcedata(1,1:end);
    
    
    % 14 columns: (sub, beta_amount, beta_prob, beta_ev, ev_bin_choice1-10)
    data_mat = zeros(length(sourcedata),14);
    
    sublist = zeros(length(sourcedata),1);
    for i = 1:length(sourcedata)
        
        % put data into table
        fname = sourcedata{i};
        T = readtable(fullfile(sourcedatadir,fname),'TreatAsEmpty','NULL');
        
        % strip out irrelevant information and missed trials
        T = T(:,{'Amount','Choice','Completed','Probability'});
        goodtrials = T.Choice < 2 & ~isnan(T.Choice);
        T = T(goodtrials,:);
        T.zAmount = zscore(T.Amount);
        T.zProbability = zscore(T.Probability);
        
        dsa = T;
        modelspec = 'Choice ~ zAmount + zProbability';
        mdl = fitglm(dsa,modelspec,'Distribution','binomial');
                
        % add expected value to the table and bin values
        T.ev = T.Amount .* T.Probability;
        bin1 = T.ev < prctile(T.ev,10);
        bin2 = T.ev >= prctile(T.ev,10) & T.ev < prctile(T.ev,20);
        bin3 = T.ev >= prctile(T.ev,20) & T.ev < prctile(T.ev,30);
        bin4 = T.ev >= prctile(T.ev,30) & T.ev < prctile(T.ev,40);
        bin5 = T.ev >= prctile(T.ev,40) & T.ev < prctile(T.ev,50);
        bin6 = T.ev >= prctile(T.ev,50) & T.ev < prctile(T.ev,60);
        bin7 = T.ev >= prctile(T.ev,60) & T.ev < prctile(T.ev,70);
        bin8 = T.ev >= prctile(T.ev,70) & T.ev < prctile(T.ev,80);
        bin9 = T.ev >= prctile(T.ev,80) & T.ev < prctile(T.ev,90);
        bin10 = T.ev >= prctile(T.ev,90);
        T.ev_binned = bin1 + bin2*2 + bin3*3 + bin4*4 + bin5*5 + bin6*6 + bin7*7 + bin8*8 + bin9*9 + bin10*10;
        
        % extract subject number from file name
        subnum = str2double(fname(3:6));
        
        
        % 14 columns: (sub, beta_amount, beta_prob, beta_ev, ev_bin_choice1-10)
        data_mat(i,1) = subnum;
        data_mat(i,2) = mdl.Coefficients.Estimate(2);
        data_mat(i,3) = mdl.Coefficients.Estimate(3);
        %data_mat(i,4) = mdl.Coefficients.Estimate(4);
        data_mat(i,4) = 0;
        for b = 1:10
            data_mat(i,b+4) = mean(T.Choice(T.ev_binned == b));
        end
    end
    
    % plot betas for logistic regression
    betas_mean = mean(data_mat(:,2:4));
    betas_se = std(data_mat(:,2:4))/sqrt(length(data_mat));
    x = 1:3;
    figure1 = figure('Name',['Logistic Regression: ' domain]);
    axes1 = axes('Parent',figure1);
    hold(axes1,'on');
    bar(x,betas_mean)
    hold on
    er = errorbar(x,betas_mean,betas_se,betas_se);
    er.Color = [0 0 0];
    er.LineStyle = 'none';
    hold off
    xlabel('Condition');
    title(['Logisitic Regression: ' domain ]);
    ylabel('Beta Coef');
    set(axes1,'XTick',[1 2 3],'XTickLabel',...
        {'Amount','Probability','Expected Value'});
    
    
    % plot mean choice values
    choice_mean = mean(data_mat(:,5:14));
    choice_se = std(data_mat(:,5:14))/sqrt(length(data_mat));
    x = 1:10;
    figure1 = figure('Name',['Choice Data: ' domain]);
    axes1 = axes('Parent',figure1);
    hold(axes1,'on');
    hold on
    er = errorbar(x,choice_mean,choice_se,choice_se);
    er.Color = [0 0 0];
    %hline(.5)
    ylim(axes1,[0 1]);
    hold off
    xlabel('Expected Value Bin');
    title(['Effort as function of Expected Value: ' domain ]);
    ylabel('Prop. Accept Hard');
    
    
end


