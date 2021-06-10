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
domain = 0;

% grab source data
% JW note [2021-06-09]: data currently on my personal computer, need to
% move it to the Smith lab computer
if domain == 0
    sourcedatadir = '/Users/jameswyngaarden/Desktop/EEfRT_SocialVsNonSocial/monetary/';
else
    if domain == 1
    sourcedatadir = '/Users/jameswyngaarden/Desktop/EEfRT_SocialVsNonSocial/social/';
    end
end
sourcedata = dir([sourcedatadir '*.xls*']);
sourcedata = struct2cell(sourcedata);
sourcedata = sourcedata(1,1:end);

% %test file for hBayesDM
% fid2 = fopen('indata_hBayesDM.tsv','w');
% fprintf(fid2,'subjID\tchoice\toutcome\n');
% datahb = zeros(100,3); % make a matrix with arbitrary number of rows

sublist = zeros(length(sourcedata),1); % add subs here, print later
idx = 0;
for i = 1:length(sourcedata)
    
    fname = sourcedata{i};
    data = readmatrix(fullfile(sourcedatadir,fname));

    % get file name parts
    fname_split = split(fname,'_');
    condition_str = fname_split{1};
    subnum_str = fname_split{2};
    
    % remove .xlsx from file name
    fname_split2 = split(subnum_str,'.');
    subnum_str = fname_split2{1};
    
    % make subject list for later
    sublist(i,1) = str2double(subnum_str);
    
    % build output file with header
    %outname = [subnum_str '_' condition_str '.csv'];
    %outfile = fullfile(maindir,'data',outname);
    %cHeader = 'subject,trial,response,reward,trial_type,accuracy';
    %fid = fopen(outfile,'w');
    %fprintf(fid,'%s\n',cHeader);
    %fclose(fid);
    
     % write data to end of file
     %dlmwrite(outfile,data,'-append');
     
%     % convert data for hBayesDM (this isn't correct, but testing/learning)
%     for t = 1:length(data)
%         if data(t,3) == -99 || data(t,4) == 0 % skip misses and neutral outcomes
%             continue
%         end
%         idx = idx + 1;
%         if strcmp(condition_str,'money')
%             datahb(idx,1) = str2double(subnum_str) + 100; %differentiate file name
%         else
%             datahb(idx,1) = str2double(subnum_str);
%         end
%         % choice
%         c = data(t,3);
%         if data(t,5) == 1 % choosing between 2 and 3 (negative)
%             if c == 2
%                 datahb(idx,2) = 1; % bad
%             elseif c == 3
%                 datahb(idx,2) = 2; % good
%             end
%         elseif data(t,5) == 0 % choosing between 1 and 3 (positive)
%             if c == 1
%                 datahb(idx,2) = 2; % good
%             elseif c == 3
%                 datahb(idx,2) = 1; % bad
%             end
%         end
%         datahb(idx,3) = data(t,4);
%     end
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
% dlmwrite('indata_hBayesDM.tsv',datahb,'delimiter','\t','-append');