% PRE-PROCESSING SCRIPT FOR BEST-CHOICE FACE TASK

% VERSION 1 (created July 2022)

% Part of the Optimal Stopping Problems Project

%% IMPORTANT NOTES %%

% I save 3 different types of mat files. 
% First type: mat files with phase 1 data (trial, block, item, rate, rt)
% Second type: mat files with phase 2 block data (trial, block, number of samples, chosen face, rank (based on rate))
% Third type: mat files with phase 2 sequence data (trial, block, every sample, every item/price, rts for each sequence)

% TOTAL MAT FILES for each subject: 62
% PHASE 1 MAT FILES: 20
% PHASE 2 BLOCK MAT FILES: 2
% PHASE 2 SEQUENCE MAT FILES: 40

% TODO:
% 1. Adapt the ideal observer model to run it through this script 
% 2. adapt model fitting to run it through this script 

%% INIT LOAD DATA %%

% get paths and define vars
startpath               = '/Users/christinadelta/githubstuff/rhul_stuff/SocialDeMa/';
resultspath             = fullfile(startpath, 'experiments', 'results');
task                    = 'face';
subpath                 = fullfile(resultspath, task);
phase                   = 1; % change this when extracting phase 2 data
session                 = 1; % change this when extracting phase 2 data

subs                    = dir(fullfile(resultspath, task, '*sub*'));
% nsubs                   = length(subs);
nsubs                   = 1;
phase1_blocks           = 20;
% phase1_blocktrials      = 40;
phase1_totaltrials      = phase1_blocks*phase1_blocktrials;

phase2_totaltrials      = 40; 
phase2_blocktrials      = 20;
phase2_blocks           = 2;
respoptions             = 2; % accept vs decline
thisequence             = nan(1,phase2_blocktrials); % every sequence has different number of draws
temp                    = 0;

% only keep subnames (not sure if this will be used)
subname                 = {subs.name};

%% EXTRACT AND SAVE PHASE 1 DATA %%

% loop over subs
for subI = 1:nsubs
    
    fprintf('loading economic phase 1 data\n')  
    subject = subs(subI).name;
    subdir  = fullfile(resultspath, task,subject);
    fprintf('\t reading data from subject %d\n',subI); 
    
    for blockI = 1:phase1_blocks
        
        fprintf('\t\t loading block %d\n\n',blockI);
        subFile = fullfile(subdir, sprintf('subject_%02d_task_%s_block_%02d_ses_%02d_phase_%02d_logs.mat',subI, task, blockI,session,phase));
        load(subFile)
        
        % extract number of this block trials
        phase1_blocktrials = size(logs.trials,2);
        
        for trial = 1:phase1_blocktrials
            
            indx = ((blockI -1)*phase1_blocktrials) + trial;  
            
            % extract info
            subj(indx)      = subI;
            trialno(indx)   = logs.trials(trial).trialNb;
            blockno(indx)   = logs.trials(trial).block;
            thisitem(indx)  = logs.trials(trial).thisitem;
            rate(indx)      = logs.trials(trial).response;
            rt(indx)        = logs.trials(trial).rt;
            
        end % end of trials loop

    end % end of block loop
     
end % end of subjects loop

% add data in one matrix
rating_data = [subj' blockno' trialno' thisitem' rate' rt'];

% % save matrix in csv format for r and python
% csvwrite('face_phase1_data.csv', rating_data)

clear subj trialno blockno thisitem rate rt session phase indx

%%  GET AVERAGE OF EACH RATING %% 

% Here for each participant we will store in s cell mean ratings
% loop over subjects 
for sub = 1:nsubs
    
    tmpsub                  = find(rating_data(:,1) == sub);
    sub_ph1_data            = rating_data((tmpsub),:);
    
    % how many unique items/faces where there to rate?
    uitems                  = length(unique(sub_ph1_data(:,4)));
    
    % init averaged ratings array
    subrate                 = nan(uitems,1);
    subitems                = [1:uitems]'; % array [1:400]
    
    for i = 1:uitems
        
        tmpitem             = find(sub_ph1_data(:,4) == i);
        tmprate             = sub_ph1_data((tmpitem),5);
        
        % average rates for item i
        subrate(i,1)        = mean(tmprate);
        subrate(i,2)        = i;
        
    end % end of unique items loop
    
    allsubs_ratings{1,sub}  = subrate;
    
    clear subrate tmpitem tmprate tmpsub subratings subitems
 
end % end of subjects loop




