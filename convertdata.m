% Extract data from xlsx and xls files and convert to csv files
%
% Also output one contenated tsv file for hbayesDM testing
% NB: hbayesDM model needs a different task structure with just two options
% so remember this is just for testing/learning purposes.
%
% 2019-12-21: created by DVS (david.v.smith@temple.edu)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 2021-06-09: updated by JBW3 (james.wyngaarden@temple.edu)
%
% First goals for update: 
% a.) loop through all of your participants
% b.) read in the data files

%% Script

clear;
maindir = pwd;

domains = {'monetary', 'social'};

m_length = 0;

for d = 1:length(domains)
    domain = domains{d};
    
    % build path for data and create list of files
    sourcedatadir = fullfile(maindir,'data',domain);
    sourcedata = dir([sourcedatadir '/*.xls*']);
    sourcedata = struct2cell(sourcedata);
    sourcedata = sourcedata(1,1:end);
    
    % make matrix to store data
    % 6 columns: domain, sub, % hard choices overall, prob 12, 50, 88
    if d == 1
        data_mat = zeros(length(sourcedata),6);
    end
       
    sublist = zeros(length(sourcedata),1); % add subs here, print later
    idx = 0;
    for i = 1:length(sourcedata)
        
        % extract data
        fname = sourcedata{i};
        data = readmatrix(fullfile(sourcedatadir,fname));
        reward = data(:,21);
        choice = data(:,22);
        reward_probability = data(:,28);
        expected_value = data(:,21) .* data(:,28) ./ 100;
        effort_data = [reward choice reward_probability expected_value];
        
        % remove misses
        effort_data(effort_data(:,2)==2,:) = [];
        
        % calculate % hard overall
        percent_hard_overall = nansum(effort_data(:,2))/(length(effort_data(:,2))-sum(isnan(effort_data(:,2))));
        
        % separate probabilities
        prob_12 = effort_data;
        prob_12(prob_12(:,3)~=12,:) = [];
        
        prob_50 = effort_data;
        prob_50(prob_50(:,3)~=50,:) = [];
        
        prob_88 = effort_data;
        prob_88(prob_88(:,3)~=88,:) = [];
        
        % percentage hard for each probability
        percent_hard_prob12 = sum(prob_12(:,2))/(length(prob_12(:,1)));
        percent_hard_prob50 = sum(prob_50(:,2))/(length(prob_50(:,1)));
        percent_hard_prob88 = sum(prob_88(:,2))/(length(prob_88(:,1)));
        
        % average expected value
        expected_value_avg = mean(expected_value);
        
        % get file name parts
        fname_split = split(fname,'_');
        condition_str = fname_split{1};
        subnum_str = fname_split{2};
        
        % remove .xlsx from file name
        fname_split2 = split(subnum_str,'.');
        subnum_str = fname_split2{1};
        
        % fill in data_mat
        if d == 1
            data_mat(i,1) = d;
            data_mat(i,2) = str2double(subnum_str);
            data_mat(i,3) = percent_hard_overall;
            data_mat(i,4) = percent_hard_prob12;
            data_mat(i,5) = percent_hard_prob50;
            data_mat(i,6) = percent_hard_prob88;
            data_mat(i,7) = expected_value_avg; % something wrong with this (NaNs)
            m_length = m_length + 1;
        else 
            if d == 2
                data_mat((i+m_length),1) = d;
                data_mat((i+m_length),2) = str2double(subnum_str);
                data_mat((i+m_length),3) = percent_hard_overall;
                data_mat((i+m_length),4) = percent_hard_prob12;
                data_mat((i+m_length),5) = percent_hard_prob50;
                data_mat((i+m_length),6) = percent_hard_prob88;
                data_mat((i+m_length),7) = expected_value_avg; % something wrong with this (NaNs)
            end
        end
        
        %keyboard
        
        % build output file with header
        %outname = [subnum_str '_' condition_str '.csv'];
        %outfile = fullfile(maindir,'data',outname);
        %cHeader = 'subject,trial,response,reward,trial_type,accuracy';
        %fid = fopen(outfile,'w');
        %fprintf(fid,'%s\n',cHeader);
        %fclose(fid);
        
        disp(domain);
        disp(subnum_str);
    end
    
end

% separate monetary and social domains from data_mat into their own
% matrices
data_mat_m = data_mat(data_mat(:,1)==1,:);
data_mat_s = data_mat(data_mat(:,1)==2,:);

% resize monetary so to remove key and sub with no social data, then
% concatenate monetary and social matrices
data_mat_m2 = data_mat_m
data_mat_m2(4,:) = []; %removes sub 1004
data_mat_m2(19,:) = []; %removes "key"
data_mat2 = [data_mat_m2 data_mat_s];


keyboard 

%% t-test: does proportion of hard-task choices overall differ between monetary and social domains?
[h,p,ci,stats] = ttest(data_mat2(:,3),data_mat2(:,10));
disp(h);
disp(p);
disp(ci);
disp(stats);

% bar graph
domain_hard_avgs = [];
domain_hard_avgs(1,1:2) = [mean(data_mat2(:,3)), mean(data_mat2(:,10))];
bar(domain_hard_avgs)
title('Proportion of hard-task choices')
xlabel('Monetary (1) Social (2)');

%% anova: does proportion of hard-task choices differ by reward probability?
% monetary
money_prob = data_mat_m(:,4:6)
[p,tbl,stats] = anova1(money_prob);

% monetary bar graph
money_prob_avgs = [];
money_prob_avgs(1,1:3) = [mean(money_prob(:,1)), mean(money_prob(:,2)), mean(money_prob(:,3))];
figure;
bar(money_prob_avgs)
title('Proportion of hard-task choices per reward probability in the monetary domain')
xlabel('12% (1)   50% (2)   88% (3)');

% social
social_prob = data_mat_s(:,4:6)
[p,tbl,stats] = anova1(social_prob);

% social bar graph
social_prob_avgs = [];
social_prob_avgs(1,1:3) = [nanmean(social_prob(:,1)), nanmean(social_prob(:,2)), nanmean(social_prob(:,3))];
figure;
bar(social_prob_avgs)
title('Proportion of hard-task choices per reward probability in the social domain')
xlabel('12% (1)   50% (2)   88% (3)');

% monetary and social together
money_prob2 = data_mat_m2(:,4:6)
total_prob = [money_prob2 social_prob];
[p,tbl,stats] = anova1(total_prob);

total_prob_avgs = [];
total_prob_avgs = [money_prob_avgs social_prob_avgs];
figure;
bar(total_prob_avgs)
title('Proportion of hard-task choices per reward probability in both domains')
xlabel('M12% (1)   M50% (2)   M88% (3)   S12% (4)   S50% (5)   S88% (6)');




%%

%plot(data_mat_m(:,7), data_mat_m(:,3));
% dlmwrite('indata_hBayesDM.tsv',datahb,'delimiter','\t','-append');

