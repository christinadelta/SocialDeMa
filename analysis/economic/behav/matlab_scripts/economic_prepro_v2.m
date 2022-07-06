% PRE-PROCESSING SCRIPT FOR BEST-CHOICE ECONOMIC TASK

% Part of the Optimal Stopping Problems Project

%% IMPORTANT NOTE %%

% I save 3 different types of mat files. 
% First type: mat files with phase 1 data (trial, block, item, rate, rt)
% Second type: mat files with phase 2 block data (trial, block, number of samples, chosen price/item, rank, reward and balance)
% Third type: mat files with phase 2 sequence data (trial, block, every sample, every item/price, rts for each sequence)

% TOTAL MAT FILES for each subject: 62
% PHASE 1 MAT FILES: 20
% PHASE 2 BLOCK MAT FILES: 2
% PHASE 2 SEQUENCE MAT FILES: 40

%% INIT LOAD DATA %%

% get paths and define vars
startpath               = '/Users/christinadelta/githubstuff/rhul_stuff/SocialDeMa/';
resultspath             = fullfile(startpath, 'experiments', 'results');
task                    = 'economic';
subpath                 = fullfile(resultspath, task);
phase                   = 1; % change this when extracting phase 2 data
session                 = 1; % change this when extracting phase 2 data

subs                    = dir(fullfile(resultspath, task, '*sub*'));
nsubs                   = length(subs);

phase1_blocks           = 20;
phase1_blocktrials      = 40;
phase1_totaltrials      = phase1_blocks*phase1_blocktrials;

phase2_totaltrials      = 40; 
phase2_blocktrials      = 20;
phase2_blocks           = 2;
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
        
        for trial = 1:phase1_blocktrials
            
            indx            = ((blockI -1)*phase1_blocktrials) + trial;  
            
            subj(indx)      = subI;
            trialno(indx)   = logs.trials(trial).trialNb;
            blockno(indx)   = logs.trials(trial).block;
            thisitem(indx)  = logs.trials(trial).thisitem;
            thisprice(indx) = logs.trials(trial).thisprice;
            rate(indx)      = logs.trials(trial).response;
            rt(indx)        = logs.trials(trial).rt;
            
        end % end of trials loop
        
    end % end of blocks loop
 
end % end of subjects loop

% add data in one matrix
rating_data = [subj' blockno' trialno' thisitem' thisprice' rate' rt'];

% % save matrix in csv format for r and python
% csvwrite('economic_phase1_data.csv', phase1_data)

clear subj trialno blockno thisitem thisprice rate rt session phase indx


%%  GET AVERAGE OF EACH RATING %% 

% loop over subjects 
for sub = 1:nsubs 
    
    tmpsub                  = find(rating_data(:,1) == sub);
    subratings              = rating_data((tmpsub),:);
    
    % how many unique prices where there to rate?
    uprice                  = length(unique(subratings(:,4)));
    
    % init averaged ratings array
    subrate                 = nan(uprice,1);
    
    % loop over unique prices
    for i = 1:uprice
        
        tmpitem             = find(subratings(:,4) == i);
        tmprate             = subratings((tmpitem),6);
        
        % average prices for item i and store
        subrate(i,1)        = mean(tmprate);
        
    end % end of unique prices loop
    
    % save this subject averaged ratings in cell
    allsubs_ratings{1,sub}  = subrate;
    
    % clear subrate tmpitem tmprate uprice tmpsub subratings
end % end of subjects loop


%% EXTRACT AND SAVE PHASE 2 BLOCK DATA %%

session     = 2;
phase       = 2;

% loop over subs
for subI = 1:nsubs
    
    fprintf('loading economic phase 2 block data\n')  
    subject = subs(subI).name;
    subdir  = fullfile(resultspath, task,subject);
    fprintf('\t reading data from subject %d\n',subI); 
    
    for blockI = 1:phase2_blocks
        
        fprintf('\t\t loading block %d\n\n',blockI);
        subFile = fullfile(subdir, sprintf('subject_%02d_task_%s_block_%02d_ses_%02d_phase_%02d_blocktrials_logs.mat',subI,task,blockI,session,phase));
        load(subFile)
        
        for trial = 1:phase2_blocktrials
            
            indx                        = ((blockI -1)*phase2_blocktrials) + trial; 
            
            subj(indx)                  = subI;
            blockno(indx)               = blockI;
            trialno(indx)               = logs.blocktrials(trial).trialnumber;
            numsamples(indx)            = logs.blocktrials(trial).numsamples;
            thisitem(indx)              = logs.blocktrials(trial).chosenitem;
            thisprice(indx)             = logs.blocktrials(trial).chosenprice;
            thisrank(indx)              = logs.blocktrials(trial).rank;
            reward(indx)                = logs.blocktrials(trial).reward;
            balance(indx)               = logs.blocktrials(trial).balance;
            
            sequences{1,subI}{1,indx}   = logs.blocktrials(trial).sequence;
            
             
        end % end of trials loop

    end % end of blocks loop
    
end % end of subjects loop

% add data in one matrix
phase2_blockdata = [subj' blockno' trialno' numsamples' thisitem' thisprice' thisrank' reward' balance'];

% % save matrix in csv format for r and python
% csvwrite('economic_phase2_blockdata.csv', phase2_blockdata)

clear subj trialno blockno thisitem thisprice numsamples thisrank indx

%% EXTRACT AND SAVE PHASE 2 SEQUENCE DATA %%

% load subject
for subI = 1:nsubs
    
    fprintf('loading phase 2 sequence data\n')  
    subject     = subs(subI).name;
    subdir      = fullfile(resultspath, task,subject);
    fprintf('\t reading data from subject %d\n',subI); 
    
    % load block
    for blockI = 1:phase2_blocks
        
        fprintf('\t\t loading phase 2 block %d\n\n',blockI);
        
        % load trial
        for trial = 1:phase2_blocktrials
            
            % load sequence file
            fprintf('\t\t loading sequence %d\n\n',trial);
            subFile = fullfile(subdir, sprintf('subject_%02d_task_%s_block_%02d_trial_%02d_ses_%02d_phase_%02d_logs.mat',subI, task, blockI, trial,session, phase));
            load(subFile)                     
            
            % how many draws in that sequence?
            thisequence(trial)          = length(logs.trialsamples);
            
            for i = 1:thisequence(trial)
                
                index                   = temp + i;    
                
                subj(index)             = subI;
                blockno(index)          = blockI;
                thisample(index)        = i;
                trialno(index)          = logs.trialsamples(i).trialnumber;
                thisitem(index)         = logs.trialsamples(i).thisitem;
                thisprice(index)        = logs.trialsamples(i).thisprice;
                rt(index)               = logs.trialsamples(i).rt;
                
            end % end of sequence loop
            
            temp                        = temp + thisequence(trial); % update temp

        end % end of trials loop
    end % end of block loop
end % end of subjects loop

% add all sequence data in one matrix
phase2_sequence_data = [subj' blockno' trialno' thisample' thisitem' thisprice' rt'];

% % save matrix in csv format for r and python
csvwrite('economic_sequencedata.csv', phase2_sequence_data)

