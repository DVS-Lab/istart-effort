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
% 1020, 1021, 1023 <- social subjects on different reward scale

%% Script

clear;
scriptsdir = pwd;
cd ..
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
        
        percent_easy_overall = choice(nansum(choice(:,1)==1))/choice(nansum(choice(:,1)));
        
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
        
        %disp(domain);
        %disp(subnum_str);
    end
end
  
% separate monetary and social domains from data_mat into their own matrices
data_mat_monetary = data_mat(data_mat(:,1)==1,:);
data_mat_social = data_mat(data_mat(:,1)==2,:);

% resize monetary and concatenate monetary and social matrices
%data_mat_monetary(4,:) = []; %removes sub 1004
%data_mat_monetary = data_mat_monetary(1:15,:);
%data_mat = [data_mat_monetary(1:15,:) data_mat_social];

writematrix(data_mat, 'data_mat.csv'); 

%%%%%%%%%%%%%%%%%%% Analyses/plots begin here %%%%%%%%%%%%%%%%%%%

%% t-test: does proportion of hard-task choices overall differ between monetary and social domains?

% organize data
domain_hard_avgs = [];
domain_hard_avgs(1,1:2) = [mean(data_mat_monetary(:,3)), mean(data_mat_social(:,3))];
sem(1,1) = std(data_mat_monetary(:,3))/sqrt(length(data_mat_monetary(:,3)));
sem(1,2) = std(data_mat_social(:,3))/sqrt(length(data_mat_social(:,3)));
x = 1:2;

% bar graph
figure1 = figure('Name','Greater Effort in Monetary Domain');
axes1 = axes('Parent',figure1);
hold(axes1,'on');
bar(x(1,1),domain_hard_avgs(1,1),'FaceColor','g');
bar(x(1,2), domain_hard_avgs(1,2),'FaceColor','b');
hold on
er = errorbar(x,domain_hard_avgs,sem,sem);
er.Color = [0 0 0];
er.LineStyle = 'none';
hold off
xlabel('Domain');
title('Greater Effort in Monetary Domain');
ylabel('Proportion of hard-task choices');
set(axes1,'XTick',[1 2],'XTickLabel',...
    {'Monetary', 'Social'});

% stats
disp("t-test: does proportion of hard-task choices overall differ between monetary and social domains?");
[h,p,ci,stats] = ttest(data_mat_monetary(:,3),data_mat_social(:,3));
disp(h);
disp(p);
disp(ci);
disp(stats);

%% anova: does proportion of hard-task choices differ by reward probability?

% organize data
money_prob = data_mat_monetary(:,4:6);
money_prob_avgs = [];
prob_avgs(1,1) = mean(money_prob(:,1));
prob_avgs(1,3) = mean(money_prob(:,2));
prob_avgs(1,5) = mean(money_prob(:,3));
social_prob = data_mat_social(:,4:6);
social_prob_avgs = [];
prob_avgs(1,2) = nanmean(social_prob(:,1));
prob_avgs(1,4) = nanmean(social_prob(:,2));
prob_avgs(1,6) = nanmean(social_prob(:,3));

sem(1,1) = std(social_prob(:,1))/sqrt(length(social_prob(:,1)));
sem(1,3) = std(social_prob(:,2))/sqrt(length(social_prob(:,2)));
sem(1,5) = std(social_prob(:,3))/sqrt(length(social_prob(:,3)));
sem(1,2) = std(social_prob(:,1))/sqrt(length(social_prob(:,1)));
sem(1,4) = std(social_prob(:,2))/sqrt(length(social_prob(:,2)));
sem(1,6) = std(social_prob(:,3))/sqrt(length(social_prob(:,3)));

% bar graph
x = 1:6;
figure1 = figure('Name','Greater Impact of Probability in Monetary Domain');
axes1 = axes('Parent',figure1);
hold(axes1,'on');
bar(x(1,1), prob_avgs(1,1),'FaceColor','g');
bar(x(1,2), prob_avgs(1,2),'FaceColor','b');
bar(x(1,3), prob_avgs(1,3),'FaceColor','g');
bar(x(1,4), prob_avgs(1,4),'FaceColor','b');
bar(x(1,5), prob_avgs(1,5),'FaceColor','g');
bar(x(1,6), prob_avgs(1,6),'FaceColor','b');
hold on
er = errorbar(x,prob_avgs,sem,sem);
er.Color = [0 0 0];
er.LineStyle = 'none';
hold off
xlabel('Probability');
title('Greater Impact of Probability in Monetary Domain');
ylabel('Proportion of hard-task choices');
set(axes1,'XTick',[1 2 3 4 5 6],'XTickLabel',...
    {'12% Mon', '12% Soc', '50% Mon', '50% Soc', '88% Mon', '88% Soc'});

% stats
[~,~,~] = anova1(money_prob);
[~,~,~] = anova1(social_prob);

%% anova: does proportion of hard-task choices differ by reward magnitude?

% organize data
money_mag = data_mat_monetary(:,7:9);
mag_avgs(1,1) = mean(money_mag(:,1));
mag_avgs(1,3) = mean(money_mag(:,2));
mag_avgs(1,5) = mean(money_mag(:,3));
social_mag = data_mat_social(:,7:9);
mag_avgs(1,2) = nanmean(social_mag(:,1));
mag_avgs(1,4) = nanmean(social_mag(:,2));
mag_avgs(1,6) = nanmean(social_mag(:,3));

sem(1,1) = nanstd(social_mag(:,1))/sqrt(length(social_mag(:,1)));
sem(1,3) = nanstd(social_mag(:,2))/sqrt(length(social_mag(:,2)));
sem(1,5) = nanstd(social_mag(:,3))/sqrt(length(social_mag(:,3)));
sem(1,2) = nanstd(social_mag(:,1))/sqrt(length(social_mag(:,1)));
sem(1,4) = nanstd(social_mag(:,2))/sqrt(length(social_mag(:,2)));
sem(1,6) = nanstd(social_mag(:,3))/sqrt(length(social_mag(:,3)));

% bar graph
x = 1:6;
figure1 = figure('Name','Greater Impact of Reward Magnitude in Monetary Domain');
axes1 = axes('Parent',figure1);
hold(axes1,'on');
bar(x(1,1), mag_avgs(1,1),'FaceColor','g');
bar(x(1,2), mag_avgs(1,2),'FaceColor','b');
bar(x(1,3), mag_avgs(1,3),'FaceColor','g');
bar(x(1,4), mag_avgs(1,4),'FaceColor','b');
bar(x(1,5), mag_avgs(1,5),'FaceColor','g');
bar(x(1,6), mag_avgs(1,6),'FaceColor','b');
hold on
er = errorbar(x,mag_avgs,sem,sem);
er.Color = [0 0 0];
er.LineStyle = 'none';
hold off
xlabel('Reward Magnitude');
title('Greater Impact of Reward Magnitude in Monetary Domain');
ylabel('Proportion of hard-task choices');
set(axes1,'XTick',[1 2 3 4 5 6],'XTickLabel',...
    {'Low Mon', 'Low Soc', 'Med Mon', 'Med Soc', 'High Mon', 'High Soc'});

% stats
[~,~,~] = anova1(money_mag);
[~,tbl,stats] = anova1(social_mag);

%% Scatterplots for proportion of hard-task choices and BDI / BSMAS for each domain

% organize data
%survey_data = readmatrix(fullfile(maindir,'data/BDI&BSMAS.xlsx'));
%BDI = survey_data(:,2);
%BSMAS = survey_data(:,3);
%y1 = data_mat_monetary(:,3);
%y2 = data_mat_social(:,3);
%
%% scatterplot for BDI
%figure1 = figure('Name','Negative Relationship between Depression and Effort in Monetary Domain');
%axes1 = axes('Parent',figure1);
%hold(axes1,'on');
%p1 = polyfit(BDI,y1,1);
%p2 = polyfit(BDI,y2,1);
%px = [min(BDI) max(BDI)];
%py1 = polyval(p1,px);
%py2 = polyval(p2,px);
%scatter(BDI,y1,'filled','g');
%hold on
%plot(px,py1,'LineWidth',2,'Color','g');
%scatter(BDI,y2,'filled','b');
%plot(px,py2,'LineWidth',2,'Color','b');
%hold off
%xlabel('BDI Score');
%title('Negative Relationship between Depression and Effort in Monetary Domain');
%ylabel('Proportion of hard-task choices');
%
%% stats
%% BDI and effort: monetary domain
%disp("BDI and effort: monetary domain (r, p)");
%[r,p] = corrcoef(BDI,y1);
%disp(r(2,1));
%disp(p(2,1));
%% BDI and effort: social domain
%disp("BDI and effort: social domain (r, p)");
%[r,p] = corrcoef(BDI,y2);
%disp(r(2,1));
%disp(p(2,1));
%
%% Reorganize to remove NaNs from BSMAS
%BSMAS(:,2) = data_mat_monetary(:,3);
%BSMAS(:,3) = data_mat_social(:,3);
%BSMAS(isnan(BSMAS(:,1)),:) = [];
%y1 = BSMAS(:,2);
%y2 = BSMAS(:,3);
%
%% scatterplot for BSMAS
%figure1 = figure('Name','Negative Relationship between Social Media Addiction and Effort in Monetary Domain');
%axes1 = axes('Parent',figure1);
%hold(axes1,'on');
%p1 = polyfit(BSMAS(:,1),y1,1);
%p2 = polyfit(BSMAS(:,1),y2,1);
%px = [min(BSMAS(:,1)) max(BSMAS(:,1))];
%py1 = polyval(p1,px);
%py2 = polyval(p2,px);
%scatter(BSMAS(:,1),y1,'filled','g');
%hold on
%plot(px,py1,'LineWidth',2,'Color','g');
%scatter(BSMAS(:,1),y2,'filled','b');
%plot(px,py2,'LineWidth',2,'Color','b');
%hold off
%xlabel('BSMAS Score');
%title('Negative Relationship between Social Media Addiction and Effort in Monetary Domain');
%ylabel('Proportion of hard-task choices');
%
%% stats
%[r,p] = corrcoef(BSMAS(:,1),y1);
%disp(r(2,1));
%disp(p(2,1));
%disp("BSMAS and effort: social domain");
%[r,p] = corrcoef(BSMAS(:,1),y2);
%disp(r(2,1));
%disp(p(2,1));
%
%
%% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% %% plotting prop of hard choices for each expected value, non-binned
% 
% fontSize1 = 18;
% fontSize2 = 15;
% 
% % monetary
% ev_data_monetary(1:4,:) = []; % remove practice trials
% 
% % create df with two columns: ev & proportion of hard task choices
% ev_probs_monetary = ev_data_monetary(:,1:2);
% for l = 1:length(ev_data_monetary(:,1))
%     ev_probs_monetary(l,1)=ev_data_monetary(l,1);
%     ev_probs_monetary(l,2)=nanmean(ev_data_monetary(l,2:length(ev_data_monetary(1,:))));
% end 
% ev_probs_monetary = sortrows(ev_probs_monetary);
% 
% % plot
% figure;
% plot(ev_probs_monetary(:,1),ev_probs_monetary(:,2))
% title('Proportion of hard-task choices for each expected value in the monetary domain', 'FontSize', fontSize1)
% xlabel('Expected Value (dollars)', 'FontSize', fontSize2)
% ylabel('Proportion of hard-task choices', 'FontSize', fontSize2);
% 
% [R,P] = corrcoef(ev_probs_monetary);
% disp(R);
% disp(P);
% 
% % social
% ev_data_social(1:4,:) = []; % remove practice trials
% 
% % create df with two columns: ev & proportion of hard task choices
% ev_probs_social = ev_data_social(:,1:2);
% for l = 1:length(ev_data_social(:,1))
%     ev_probs_social(l,1)=ev_data_social(l,1);
%     ev_probs_social(l,2)=nanmean(ev_data_social(l,2:length(ev_data_social(1,:))));
% end 
% ev_probs_social = sortrows(ev_probs_social);
% 
% % plot
% figure;
% plot(ev_probs_social(:,1),ev_probs_social(:,2))
% title('Proportion of hard-task choices for each expected value in the social domain', 'FontSize', fontSize1)
% xlabel('Expected Value (minutes of social media time)', 'FontSize', fontSize2)
% ylabel('Proportion of hard-task choices', 'FontSize', fontSize2);
% 
% [R,P] = corrcoef(ev_probs_social);
% disp(R);
% disp(P);
% 
% 
% 
% 
% %% plotting prop of hard choices for each reward magnitude, non-binned
% 
% % NOTE: reward values get shown multiple times (i.e., the same value will
% % appear in column 1 of reward_probs more than once), need to find a way to
% % combine those, probably need to go back to the data file and sort before
% % doing any of the matrix stuff in the main script above
% 
% %monetary
% reward_data_monetary(1:4,:) = [];
% 
% reward_probs_monetary = reward_data_monetary(:,1:2);
% for l = 1:length(reward_data_monetary(:,1))
%     reward_probs_monetary(l,1)=reward_data_monetary(l,1);
%     reward_probs_monetary(l,2)=nanmean(reward_data_monetary(l,2:length(reward_data_monetary(1,:))));
% end 
% reward_probs_monetary = sortrows(reward_probs_monetary);
% 
% % social
% reward_data_social(1:4,:) = [];
% 
% reward_probs_social = reward_data_social(:,1:2);
% for l = 1:length(reward_data_social(:,1))
%     reward_probs_social(l,1)=reward_data_social(l,1);
%     reward_probs_social(l,2)=nanmean(reward_data_social(l,2:length(reward_data_social(1,:))));
% end 
% reward_probs_social = sortrows(reward_probs_social);
% 
% %%
% 
% %plot(data_mat_m(:,7), data_mat_m(:,3));
% % dlmwrite('indata_hBayesDM.tsv',datahb,'delimiter','\t','-append');

