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

clear;
maindir = pwd;

domains = {'monetary', 'social'};
for d = 1:length(domains)
    domain = domains{d};
    
    % build path for data and create list of files
    sourcedatadir = fullfile(maindir,'data',domain);
    sourcedata = dir([sourcedatadir '/*.xls*']);
    sourcedata = struct2cell(sourcedata);
    sourcedata = sourcedata(1,1:end);
    
    % make matrix to store data
    % 5 columns: sub, % hard choices overall, prob 12, 50, 88
    data_mat = zeros(length(sourcedata),5);
    
    
    sublist = zeros(length(sourcedata),1); % add subs here, print later
    idx = 0;
    for i = 1:length(sourcedata)
        
        % extract data
        fname = sourcedata{i};
        data = readmatrix(fullfile(sourcedatadir,fname));
        choice = data(:,22);
        reward_probability = data(:,28);
        effort_data = [choice reward_probability];
        
        % remove misses
        effort_data(effort_data(:,1)==2,:) = [];
        
        % calculate % hard overall
        percent_hard_overall = nansum(effort_data(:,1))/(length(effort_data(:,1))-sum(isnan(effort_data(:,1))));
        
        % separate probabilities
        prob_12 = effort_data;
        prob_12(prob_12(:,2)~=12,:) = [];
        
        prob_50 = effort_data;
        prob_50(prob_50(:,2)~=50,:) = [];
        
        prob_88 = effort_data;
        prob_88(prob_88(:,2)~=88,:) = [];
        
        % percentage hard for each probability
        percent_hard_prob12 = sum(prob_12(:,1))/(length(prob_12(:,1)));
        percent_hard_prob50 = sum(prob_50(:,1))/(length(prob_50(:,1)));
        percent_hard_prob88 = sum(prob_88(:,1))/(length(prob_88(:,1)));
        
        % get file name parts
        fname_split = split(fname,'_');
        condition_str = fname_split{1};
        subnum_str = fname_split{2};
        
        % remove .xlsx from file name
        fname_split2 = split(subnum_str,'.');
        subnum_str = fname_split2{1};
        
        % fill in data_mat
        data_mat(i,1) = d;
        data_mat(i,2) = str2double(subnum_str);
        data_mat(i,3) = percent_hard_overall;
        data_mat(i,4) = percent_hard_prob12;
        data_mat(i,5) = percent_hard_prob50;
        data_mat(i,6) = percent_hard_prob88;
        
        %keyboard
        % make subject list for later
        sublist(i,1) = str2double(subnum_str);
        
        
        
        % build output file with header
        %outname = [subnum_str '_' condition_str '.csv'];
        %outfile = fullfile(maindir,'data',outname);
        %cHeader = 'subject,trial,response,reward,trial_type,accuracy';
        %fid = fopen(outfile,'w');
        %fprintf(fid,'%s\n',cHeader);
        %fclose(fid);
        
        disp(subnum_str);
    end
    
    % % write out a sublist for later
    if domain == 0
        dlmwrite('M_sublist.txt',unique(sublist));
    else
        if domain == 1
            dlmwrite('S_sublist.txt',unique(sublist));
        end
    end
end
% dlmwrite('indata_hBayesDM.tsv',datahb,'delimiter','\t','-append');