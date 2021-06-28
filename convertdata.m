% Extract data from xlsx and xls files and convert to csv files
%
% Also output one contenated tsv file for hbayesDM testing
% NB: hbayesDM model needs a different task structure with just two options
% so remember this is just for testing/learning purposes.
%
% 2019-12-21: created by DVS (david.v.smith@temple.edu)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 2021-06-09: updated for ISTART effort task data by JBW3 (james.wyngaarden@temple.edu)

%% Script

clear;
maindir = pwd;

domains = {'monetary', 'social'};

% preallocate vars that otherwise change sizes in loop below:
m_length = 0;
reward_data_monetary = NaN(94,20);
reward_data_social = NaN(94,20);
ev_data_monetary = reward_data_monetary;
ev_data_social = reward_data_social;

% first loop through monetary domain, then social
for d = 1:length(domains)
    domain = domains{d};
    
    % build path for data and create list of files
    sourcedatadir = fullfile(maindir,'data',domain);
    sourcedata = dir([sourcedatadir '/*.xls*']);
    sourcedata = struct2cell(sourcedata);
    sourcedata = sourcedata(1,1:end);
    
    % make matrix to store data
    % 18 columns: (domain, sub, % hard choices overall, prob 12, 50, 88, prob
    % low reward, prob med reward, prob hi reward)*2 domains 
    if d == 1
        data_mat = zeros(length(sourcedata),9);
    end
       
    sublist = zeros(length(sourcedata),1);
    idx = 0;
    for i = 1:length(sourcedata)
        
        % extract data
        fname = sourcedata{i};
        data = readmatrix(fullfile(sourcedatadir,fname));
        reward = data(:,21);
        choice = data(:,22);
        reward_probability = data(:,28);
        expected_value = reward*0; % preallocating length with zeros
        for j = 1:length(data(:,21))
            expected_value(j,1) = data(j,21) .* data(j,28) ./ 100;
        end
        reward_bins = reward*0;
        effort_data = [reward reward_bins reward_probability choice expected_value];
        
        % remove NaNs & misses
        effort_data(any(isnan(effort_data),1),:) = [];
        effort_data(effort_data(:,4)==2,:) = [];
        
        % calculate % hard overall
        percent_hard_overall = nansum(effort_data(:,4))/(length(effort_data(:,4))-sum(isnan(effort_data(:,4))));
        
        % separate probabilities & find percent hard for each
        prob_12 = effort_data;
        prob_12(prob_12(:,3)~=12,:) = [];
        prob_50 = effort_data;
        prob_50(prob_50(:,3)~=50,:) = [];
        prob_88 = effort_data;
        prob_88(prob_88(:,3)~=88,:) = [];
        percent_hard_prob12 = sum(prob_12(:,4))/(length(prob_12(:,4)));
        percent_hard_prob50 = sum(prob_50(:,4))/(length(prob_50(:,4)));
        percent_hard_prob88 = sum(prob_88(:,4))/(length(prob_88(:,4)));
        
        % separate reward magnitudes into bins
        % monetary: low = 1.24-2.23; med = 2.24-3.22; high = 3.23-4.21
        if d == 1
            for f = 1:length(effort_data(:,1))
                if (1.23<(effort_data(f,1))) && ((effort_data(f,1))<2.24)
                    effort_data(f,2) = 1;
                else
                    if (2.23<(effort_data(f,1))) && ((effort_data(f,1))<3.23)
                        effort_data(f,2) = 2;
                    else
                        if (3.22<(effort_data(f,1))) && ((effort_data(f,1))<4.22)
                        effort_data(f,2) = 3;
                        end
                    end
                end
            end
        else
            
        % social: low = 8.65-15.57; med = 15.58-22.46; high = 22.47-29.37
        % note: 1020, 1021, and 1023 have different amounts for social than the
        % other participants, have been removed from data dir
            if d == 2
                for f = 1:length(effort_data(:,1))
                    if (8.64<(effort_data(f,1))) && ((effort_data(f,1))<15.58)
                        effort_data((f+m_length),2) = 1;
                    else
                        if (15.57<(effort_data(f,1))) && ((effort_data(f,1))<22.47)
                            effort_data((f+m_length),2) = 2;
                        else
                            if (22.46<(effort_data(f,1))) && ((effort_data(f,1))<29.38)
                                effort_data((f+m_length),2) = 3;
                            end
                        end
                    end
                end
            end
        end    
        
        % separate reward bins & find percent hard for each
        reward_low = effort_data;
        reward_low(reward_low(:,2)~=1,:) = [];
        reward_mid = effort_data;
        reward_mid(reward_mid(:,2)~=2,:) = [];
        reward_hi = effort_data;
        reward_hi(reward_hi(:,2)~=3,:) = [];
        percent_hard_reward_low = sum(reward_low(:,4))/(length(reward_low(:,4)));
        percent_hard_reward_mid = sum(reward_mid(:,4))/(length(reward_mid(:,4)));
        percent_hard_reward_hi = sum(reward_hi(:,4))/(length(reward_hi(:,4)));
        
        % prob of picking hard for each individual reward value:
        if d == 1
            if length(reward(:,1))>length(reward_data_monetary(:,1))
                diff = length(reward(:,1)) - length(reward_data_monetary(:,1));
                reward_data_monetary(1:length(reward(:,1)),1) = reward;
                reward_data_monetary((length(reward(:,1))-diff):(length(reward(:,1))),(2):(length(reward_data_monetary(1,:)))) = NaN;
            end
            for h = 1:length(choice(:,1)) 
            reward_data_monetary(h,1+i) = choice(h,1);
            end
        else
            if d == 2
                if length(reward(:,1))>length(reward_data_social(:,1))
                    diff = length(reward(:,1)) - length(reward_data_social(:,1));
                    reward_data_social(1:length(reward(:,1)),1) = reward;
                    reward_data_social((length(reward(:,1))-diff):(length(reward(:,1))),(2):(length(reward_data_social(1,:)))) = NaN;
                end
                for h = 1:length(choice(:,1)) 
                reward_data_social(h,1+i) = choice(h,1);
                end
            end
        end
        
        % prob of picking hard for each expected value
        if d == 1
            if length(expected_value(:,1))>length(ev_data_monetary(:,1))
                diff = length(expected_value(:,1)) - length(ev_data_monetary(:,1));
                ev_data_monetary(1:length(expected_value(:,1)),1) = expected_value;
                ev_data_monetary((length(expected_value(:,1))-diff):(length(expected_value(:,1))),(2):(length(ev_data_monetary(1,:)))) = NaN;
            end
            for h = 1:length(choice(:,1)) 
            ev_data_monetary(h,1+i) = choice(h,1);
            end
        else
            if d == 2
                if length(expected_value(:,1))>length(ev_data_social(:,1))
                    diff = length(expected_value(:,1)) - length(ev_data_social(:,1));
                    ev_data_social(1:length(expected_value(:,1)),1) = expected_value;
                    ev_data_social((length(expected_value(:,1))-diff):(length(expected_value(:,1))),(2):(length(ev_data_social(1,:)))) = NaN;
                end
                for h = 1:length(choice(:,1)) 
                ev_data_social(h,1+i) = choice(h,1);
                end
            end
        end
        
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
            data_mat(i,7) = percent_hard_reward_low;
            data_mat(i,8) = percent_hard_reward_mid;
            data_mat(i,9) = percent_hard_reward_hi;
            m_length = m_length + 1;
        else 
            if d == 2
                data_mat((i+m_length),1) = d;
                data_mat((i+m_length),2) = str2double(subnum_str);
                data_mat((i+m_length),3) = percent_hard_overall;
                data_mat((i+m_length),4) = percent_hard_prob12;
                data_mat((i+m_length),5) = percent_hard_prob50;
                data_mat((i+m_length),6) = percent_hard_prob88;
                data_mat((i+m_length),7) = percent_hard_reward_low;
                data_mat((i+m_length),8) = percent_hard_reward_mid;
                data_mat((i+m_length),9) = percent_hard_reward_hi;
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
        
        % print domain & subject ID while running
        disp(domain);
        disp(subnum_str);
    end
end
  
% separate monetary and social domains from data_mat into their own matrices
data_mat_monetary = data_mat(data_mat(:,1)==1,:);
data_mat_social = data_mat(data_mat(:,1)==2,:);

% resize monetary and concatenate monetary and social matrices
data_mat_monetary(4,:) = []; %removes sub 1004
data_mat = [data_mat_monetary(1:15,:) data_mat_social];

writematrix(data_mat, 'data_mat.csv');

keyboard 

%%%%%%%%%%%%%%%%%%% Analyses/plots begin here %%%%%%%%%%%%%%%%%%%

%% t-test: does proportion of hard-task choices overall differ between monetary and social domains?
[h,p,ci,stats] = ttest(data_mat2(:,3),data_mat2(:,12));
disp(h);
disp(p);
disp(ci);
disp(stats);

% bar graph
domain_hard_avgs = [];
domain_hard_avgs(1,1:2) = [mean(data_mat2(:,3)), mean(data_mat2(:,12))];

% standard error
sem1 = std(data_mat2(:,3))/sqrt(length(data_mat2(:,3)));
sem2 = std(data_mat2(:,12))/sqrt(length(data_mat2(:,12)));
err = [sem1 sem2];

bar(domain_hard_avgs)
title('Proportion of hard-task choices')
ylabel('proportion of hard-task choices')
xlabel('Monetary (1) Social (2)')

hold on

er = errorbar(domain_hard_avgs,err);
er.Color = [0 0 0];
er.LineStyle = 'none';

hold off

%% anova: does proportion of hard-task choices differ by reward probability?
% monetary
money_prob = data_mat_m(:,4:6);
[~,~,~] = anova1(money_prob);

% monetary bar graph
money_prob_avgs = [];
money_prob_avgs(1,1:3) = [mean(money_prob(:,1)), mean(money_prob(:,2)), mean(money_prob(:,3))];
figure;
bar(money_prob_avgs)
title('Proportion of hard-task choices per reward probability in the monetary domain')
ylabel('proportion of hard-task choices')
xlabel('12% (1)   50% (2)   88% (3)')
ylim([0 .7]);

% social
social_prob = data_mat_s(:,4:6);
[~,~,~] = anova1(social_prob);

% social bar graph
social_prob_avgs = [];
social_prob_avgs(1,1:3) = [nanmean(social_prob(:,1)), nanmean(social_prob(:,2)), nanmean(social_prob(:,3))];
figure;
bar(social_prob_avgs)
title('Proportion of hard-task choices per reward probability in the social domain')
ylabel('proportion of hard-task choices')
xlabel('12% (1)   50% (2)   88% (3)')
ylim([0 .7]);

% monetary and social together
money_prob2 = data_mat_m2(:,4:6);
total_prob = [money_prob2 social_prob];
[~,~,~] = anova1(total_prob);

total_prob_avgs = [money_prob_avgs social_prob_avgs];
figure;
bar(total_prob_avgs)
title('Proportion of hard-task choices per reward probability in both domains')
ylabel('proportion of hard-task choices')
xlabel('M12% (1)   M50% (2)   M88% (3)   S12% (4)   S50% (5)   S88% (6)');

%% anova: does proportion of hard-task choices differ by reward magnitude?
% monetary
money_mag = data_mat_m(:,7:9);
[~,~,~] = anova1(money_mag);

% monetary bar graph
money_mag_avgs = [];
money_mag_avgs(1,1:3) = [mean(money_mag(:,1)), mean(money_mag(:,2)), mean(money_mag(:,3))];

% standard error
sem1 = std(money_mag(:,1))/sqrt(length(money_mag(:,1)));
sem2 = std(money_mag(:,2))/sqrt(length(money_mag(:,2)));
sem3 = std(money_mag(:,3))/sqrt(length(money_mag(:,3)));
err = [sem1 sem2 sem3];

figure;
bar(money_mag_avgs)
title('Proportion of hard-task choices per reward magnitude in the monetary domain')
ylabel('proportion of hard-task choices')
xlabel('low reward (1)   mid reward (2)   high reward (3)')
ylim([0 .7])

hold on

er = errorbar(money_mag_avgs,err);
er.Color = [0 0 0];
er.LineStyle = 'none';

hold off

% social
social_mag = data_mat_s(:,7:9);
social_mag(any(isnan(social_mag),1),:) = [];
[p,tbl,stats] = anova1(social_mag);

% social bar graph
social_mag_avgs = [];
social_mag_avgs(1,1:3) = [nanmean(social_mag(:,1)), nanmean(social_mag(:,2)), nanmean(social_mag(:,3))];

% standard error
sem1 = nanstd(social_mag(:,1))/sqrt(length(social_mag(:,1)));
sem2 = nanstd(social_mag(:,2))/sqrt(length(social_mag(:,2)));
sem3 = nanstd(social_mag(:,3))/sqrt(length(social_mag(:,3)));
err = [sem1 sem2 sem3];

figure;
bar(social_mag_avgs)
title('Proportion of hard-task choices per reward magnitude in the social domain')
ylabel('proportion of hard-task choices')
xlabel('low reward (1)   mid reward (2)   high reward (3)')
ylim([0 .7])

hold on

er = errorbar(social_mag_avgs,err);
er.Color = [0 0 0];
er.LineStyle = 'none';

hold off

%% plotting prop of hard choices for each expected value, non-binned

fontSize1 = 18;
fontSize2 = 15;

% monetary
ev_data_monetary(1:4,:) = []; % remove practice trials

% create df with two columns: ev & proportion of hard task choices
ev_probs_monetary = ev_data_monetary(:,1:2);
for l = 1:length(ev_data_monetary(:,1))
    ev_probs_monetary(l,1)=ev_data_monetary(l,1);
    ev_probs_monetary(l,2)=nanmean(ev_data_monetary(l,2:length(ev_data_monetary(1,:))));
end 
ev_probs_monetary = sortrows(ev_probs_monetary);

% plot
figure;
plot(ev_probs_monetary(:,1),ev_probs_monetary(:,2))
title('Proportion of hard-task choices for each expected value in the monetary domain', 'FontSize', fontSize1)
xlabel('Expected Value (dollars)', 'FontSize', fontSize2)
ylabel('Proportion of hard-task choices', 'FontSize', fontSize2);

% social
ev_data_social(1:4,:) = []; % remove practice trials

% create df with two columns: ev & proportion of hard task choices
ev_probs_social = ev_data_social(:,1:2);
for l = 1:length(ev_data_social(:,1))
    ev_probs_social(l,1)=ev_data_social(l,1);
    ev_probs_social(l,2)=nanmean(ev_data_social(l,2:length(ev_data_social(1,:))));
end 
ev_probs_social = sortrows(ev_probs_social);

% plot
figure;
plot(ev_probs_social(:,1),ev_probs_social(:,2))
title('Proportion of hard-task choices for each expected value in the social domain', 'FontSize', fontSize1)
xlabel('Expected Value (minutes of social media time)', 'FontSize', fontSize2)
ylabel('Proportion of hard-task choices', 'FontSize', fontSize2);

%% plotting prop of hard choices for each reward magnitude, non-binned

% NOTE: reward values get shown multiple times (i.e., the same value will
% appear in column 1 of reward_probs more than once), need to find a way to
% combine those, probably need to go back to the data file and sort before
% doing any of the matrix stuff in the main script above

%monetary
reward_data_monetary(1:4,:) = [];

reward_probs_monetary = reward_data_monetary(:,1:2);
for l = 1:length(reward_data_monetary(:,1))
    reward_probs_monetary(l,1)=reward_data_monetary(l,1);
    reward_probs_monetary(l,2)=nanmean(reward_data_monetary(l,2:length(reward_data_monetary(1,:))));
end 
reward_probs_monetary = sortrows(reward_probs_monetary);

% social
reward_data_social(1:4,:) = [];

reward_probs_social = reward_data_social(:,1:2);
for l = 1:length(reward_data_social(:,1))
    reward_probs_social(l,1)=reward_data_social(l,1);
    reward_probs_social(l,2)=nanmean(reward_data_social(l,2:length(reward_data_social(1,:))));
end 
reward_probs_social = sortrows(reward_probs_social);

%%

%plot(data_mat_m(:,7), data_mat_m(:,3));
% dlmwrite('indata_hBayesDM.tsv',datahb,'delimiter','\t','-append');

