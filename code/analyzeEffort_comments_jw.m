% Analyze sourcedata (.xlsx) from Effort Task. Goal is to quantify how
% willingness to exert effort in social and monetary conditions is
% influenced by value.

% Notes:
% domain is the same thing as "condition"
% 1020, 1021, 1023 <- social subjects on different reward scale.
% need one subject list


clear; close all; %clear: removes all vars from the workspace, removes functions from memory; close: closes figure windows
maindir = pwd; %pwd: present working directory (i.e., .../istart-effort/); maindir is a character var
warning off all %warning: disable or enable warning messages


% first loop through monetary domain, then social
domains = {'monetary', 'social'}; %domains: character array
for d = 1:length(domains) %d: also character array?
    domain = domains{d}; %domain: character var
    
    % build path for data and create list of files
    sourcedatadir = fullfile(maindir,'data',domain); %fullfile: builds full file name from parts; parts can be strings, character vectors, or cell arrays of character vectors
    sourcedata = dir([sourcedatadir '/*.xls*']); %dir: lists a folder; dir NAME lists the files within a folder. asterisks = wildcards
    sourcedata = struct2cell(sourcedata); %struct2cell: converst structure array to cell array
    sourcedata = sourcedata(1,1:end); %trims sourcedata down to only the first row (i.e., sub ID files); keeps all columns (i.e., all subs); sourcedata: cell array, class=char?
    
    
    % 14 columns: (sub, beta_amount, beta_prob, beta_ev, ev_bin_choice1-10)
    data_mat = zeros(length(sourcedata),14); %zeros: creates and N-by-N matrix of zeros
    
    sublist = zeros(length(sourcedata),1);
    for i = 1:length(sourcedata)
        
        % put data into table
        fname = sourcedata{i};
        % sub-1007 does not have NULL
        T = readtable(fullfile(sourcedatadir,fname),'TreatAsEmpty','NULL'); %readtable: creates a table by reading from a file; readtable(FILENAME); 'TreatAsEmpty': text which is used in the file to represent empty data, 'NULL' means empty data gets turned into NaNs
        
        % strip out irrelevant information and missed trials
        T = T(:,{'Amount','Choice','Completed','Probability'});
        
        easytrials = T.Choice == 0 & ~isnan(T.Choice);
        hardtrials = T.Choice == 1 & ~isnan(T.Choice);
        noresptrials = T.Choice == 2 & ~isnan(T.Choice);
        choicelength = ~isnan(T.Choice);
        
        easyprob = sum(easytrials)/length(choicelength);
        hardprob = sum(hardtrials)/length(choicelength);
        norespprob = sum(noresptrials)/length(choicelength);
        
        goodtrials = T.Choice < 2 & ~isnan(T.Choice) & T.Amount > 0;
        T = T(goodtrials,:);
        T.zAmount = zscore(T.Amount);
        T.zProbability = zscore(T.Probability);
        
        dsa = T;
        %modelspec = 'Choice ~ zAmount + zProbability';
        modelspec = 'Choice ~ zAmount*zProbability';
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
        
        %if subnum == 1007 && strcmp(domain,'social')
        %    keyboard
        %end
        
        % 14 columns: (sub, beta_amount, beta_prob, beta_ev, ev_bin_choice1-10)
        data_mat(i,1) = subnum;
        %data_mat(i,2) = mdl.Coefficients.Estimate(2);
        %data_mat(i,3) = mdl.Coefficients.Estimate(3);
        %data_mat(i,4) = mdl.Coefficients.Estimate(4);
        data_mat(i,2) = mdl.Coefficients.tStat(2);
        data_mat(i,3) = mdl.Coefficients.tStat(3);
        data_mat(i,4) = mdl.Coefficients.tStat(4);
        %data_mat(i,4) = 0;
        for b = 1:10
            data_mat(i,b+4) = mean(T.Choice(T.ev_binned == b));
        end
        data_mat(i,15) = easyprob;
        data_mat(i,16) = hardprob;
        data_mat(i,17) = norespprob;
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
    ylabel('t-stat');
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
    
    % compile beta stats, 4 columns: M betas se, M betas mean, S se, S mean
    if d == 1
        betas_se_both(1,1) = betas_se(1,1);
        betas_se_both(1,3) = betas_se(1,2);
        %betas_se_both(1,5) = betas_se(1,3);
        betas_mean_both(1,1) = betas_mean(1,1);
        betas_mean_both(1,3) = betas_mean(1,2);
        %betas_mean_both(1,5) = betas_mean(1,3);
    else
        betas_se_both(1,2) = betas_se(1,1);
        betas_se_both(1,4) = betas_se(1,2);
        %betas_se_both(1,6) = betas_se(1,3);
        betas_mean_both(1,2) = betas_mean(1,1);
        betas_mean_both(1,4) = betas_mean(1,2);
        %betas_mean_both(1,6) = betas_mean(1,3);
    end
    
end

% % plot both domains together
% beta_stats = beta_stats(1:2,:);
% figure, barwitherr([beta_stats(:,1);beta_stats(:,3)],[beta_stats(:,2);beta_stats(:,4)]);
% title('Logisitic Regression');
% ylabel('t-stat');
% set(axes1,'XTick',[1 2 3 4],'XTickLabel',...
%     {'Amount','Probability'});


% plot betas for logistic regression
x = 1:4;
figure1 = figure('Name','Logistic Regression: ');
axes1 = axes('Parent',figure1);
hold(axes1,'on');
bar(x,betas_mean_both)
hold on
er = errorbar(x,betas_mean_both,betas_se_both,betas_se_both);
er.Color = [0 0 0];
er.LineStyle = 'none';
hold off
xlabel('Condition');
title('Logisitic Regression: both domains');
ylabel('t-stat');
set(axes1,'XTick',[1 2 3 4],'XTickLabel',...
    {'Monetary Amount','Social Amount','Monetary Probability','Social Probability'});
% figure1(1).FaceColor = 'r';
% figure1(2).FaceColor = 'b';
% figure1(3).FaceColor = 'r';
% figure1(4).FaceColor = 'b';
% figure1(5).FaceColor = 'r';
% figure1(6).FaceColor = 'b';

