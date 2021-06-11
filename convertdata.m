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

% set domain (i.e., 0 = monetary, 1 = social)
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
        
        percent_hard_overall = sum(effort_data(:,1))/length(effort_data);
        
        keyboard
        % get file name parts
        fname_split = split(fname,'_');
        condition_str = fname_split{1};
        subnum_str = fname_split{2};
        
        % remove .xlsx from file name
        fname_split2 = split(subnum_str,'.');
        subnum_str = fname_split2{1};
        
        % make subject list for later
        sublist(i,1) = str2double(subnum_str);
        
        % feed data into matrix for hard vs. easy
        for j = 1:length(data(:,22))
            disp(j);
            if data(j:22) == 2
                data(j:22) = NaN;
            end
        end
        
        hard = nansum(data(:,22));
        
        % this works:
        % data(65,22)=NaN;
        
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