% PRE-PROCESSING SCRIPT FOR BEST-CHOICE ECONOMIC TASK

% VERSION 3 (updated July 2022)

% SECOND UPDATE: 16/01/2023

% Part of the Optimal Stopping Problems Project

%% IMPORTANT NOTES %%

% I save 3 different types of mat files. 
% First type: mat files with phase 1 data (trial, block, item, rate, rt)
% Second type: mat files with phase 2 block data (trial, block, number of samples, chosen price/item, rank, reward and balance)
% Third type: mat files with phase 2 sequence data (trial, block, every sample, every item/price, rts for each sequence)


% TOTAL MAT FILES for each subject: 62
% PHASE 1 MAT FILES: 20
% PHASE 2 BLOCK MAT FILES: 2
% PHASE 2 SEQUENCE MAT FILES: 40

% The script also deals with sequences, choiceVectors, prices, ratings to run
% the ideal observer model, also for model fitting and model comparison 

% TODO:
% 1. Adapt the ideal observer model to run it through this script 
% 2. adapt model fitting to run it through this script 

%% INIT LOAD DATA %%

% get paths and define vars
startpath               = '/Volumes/DeepSpaceStuff/optimal_stopping_data/data/';
% resultspath             = fullfile(startpath, 'experiments', 'results');
task                    = 'economic';
datatype                = 'behav'; 
subpath                 = fullfile(startpath, task, datatype);
phase                   = 1; % change this when extracting phase 2 data
session                 = 1; % change this when extracting phase 2 data

subs                    = dir(fullfile(subpath, '*sub*'));
nsubs                   = length(subs);
% nsubs                   = 5;

phase1_blocks           = 20;
phase1_blocktrials      = 40;
phase1_totaltrials      = phase1_blocks*phase1_blocktrials;

phase2_totaltrials      = 40; 
phase2_blocktrials      = 20;
phase2_blocks           = 2;
respoptions             = 2; % accept vs decline (for phase 2)
thisequence             = nan(1,phase2_blocktrials); % every sequence has different number of draws
temp                    = 0;

counter                 = 0; 
% only keep subnames (not sure if this will be used)
subname                 = {subs.name};

%% EXTRACT AND SAVE PHASE 1 DATA %%

% loop over subs
for subI = 1:nsubs
    
    fprintf('loading economic phase 1 data\n')  
    subject = subs(subI).name;
    subdir  = fullfile(subpath,subject);
    fprintf('\t reading data from subject %d\n',subI); 
    
    for blockI = 1:phase1_blocks
        
        fprintf('\t\t loading block %d\n\n',blockI);
        subFile = fullfile(subdir, sprintf('subject_%02d_task_%s_block_%02d_ses_%02d_phase_%02d_logs.mat',subI, task, blockI,session,phase));
        load(subFile)
        
        for trial = 1:phase1_blocktrials
            
            indx            = counter + ((blockI -1)*phase1_blocktrials) + trial;  
            
            subj(indx)      = subI;
            trialno(indx)   = logs.trials(trial).trialNb;
            blockno(indx)   = logs.trials(trial).block;
            thisitem(indx)  = logs.trials(trial).thisitem;
            thisprice(indx) = logs.trials(trial).thisprice;
            rate(indx)      = logs.trials(trial).response;
            rt(indx)        = logs.trials(trial).rt;
            
        end % end of trials loop
        
    end % end of blocks loop
    
    % update indx var so that it carries on after each participant
    counter                 = counter + (indx/subI); 
    
end % end of subjects loop

% add data in one matrix
rating_data                 = [subj' blockno' trialno' thisitem' thisprice' rate' rt'];

% % save matrix in csv format for r and python
% csvwrite('economic_phase1_data.csv', phase1_data)

clear subj trialno blockno thisitem thisprice rate rt session phase indx


%%  GET AVERAGE OF EACH RATING %% 

% Here for each participant we will store in s cell mean ratings
% loop over subjects 
for sub = 1:nsubs 
    
    tmpsub                  = find(rating_data(:,1) == sub);
    sub_ph1_data            = rating_data((tmpsub),:);
    
    % how many unique prices where there to rate?
    uitems                  = length(unique(sub_ph1_data(:,4)));
    
    % init averaged ratings array
    subrate                 = nan(uitems,1);
    subitems                = [1:uitems]'; % array [1:400]
    
    % loop over unique prices
    for i = 1:uitems
        
        tmpitem             = find(sub_ph1_data(:,4) == i);
        tmprate             = sub_ph1_data((tmpitem),6);
        tmpprice            = sub_ph1_data((tmpitem),5);
        
        % average prices for item i and store
        subrate(i,1)        = mean(tmprate);
        subrate(i,2)        = i;
        
        % each price is shown twise, so the tmpprice placeholder contains
        % the same price twise, just pick the first one
        subitems(i,2)       = tmpprice(1);
        
    end % end of unique prices loop
    
    % extract raw prices for each item
    
    
    % save this subject averaged ratings in cell
    allsubs_ratings{1,sub}  = subrate;
    allsubs_prices{1,sub}   = subitems;
    
    clear subrate tmpitem tmprate uprice tmpsub subratings subitems
end % end of subjects loop


%% EXTRACT AND SAVE PHASE 2 BLOCK DATA %%

session     = 2;
phase       = 2;
counter     = 0; % init counter var

% loop over subs
for subI = 1:nsubs
    
    fprintf('loading economic phase 2 block data\n')  
    subject = subs(subI).name;
    subdir  = fullfile(subpath,subject);
    fprintf('\t reading data from subject %d\n',subI); 
    
    for blockI = 1:phase2_blocks
        
        fprintf('\t\t loading block %d\n\n',blockI);
        subFile = fullfile(subdir, sprintf('subject_%02d_task_%s_block_%02d_ses_%02d_phase_%02d_blocktrials_logs.mat',subI,task,blockI,session,phase));
        load(subFile)
        
        for trial = 1:phase2_blocktrials
            
            indx                                = counter + ((blockI -1)*phase2_blocktrials) + trial; 
            id                                  = ((blockI -1)*phase2_blocktrials) + trial;
            
            subj(indx)                          = subI;
            blockno(indx)                       = blockI;
            trialno(indx)                       = logs.blocktrials(trial).trialnumber;
            numsamples(indx)                    = logs.blocktrials(trial).numsamples;
            thisitem(indx)                      = logs.blocktrials(trial).chosenitem;
            thisprice(indx)                     = logs.blocktrials(trial).chosenprice;
            thisrank(indx)                      = logs.blocktrials(trial).rank;
            reward(indx)                        = logs.blocktrials(trial).reward;
            balance(indx)                       = logs.blocktrials(trial).balance;
            
            allsubs_sequences{1,subI}{1,id}   = logs.blocktrials(trial).sequence;
            
            % create a vector of sequence responses at this point. This
            % will be based on the number of samples on every
            % trial/sequence
            t                                   = nan(numsamples(indx), respoptions); % init empty vec
            
            for s = 1:numsamples(indx) 
                
                if s < numsamples(indx) % if this is a decline response
                    t(s,1)                      = 1;
                    t(s,2)                      = 0;
                else                    % if participant accepted an option
                    t(s,1)                      = 0;
                    t(s,2)                      = 1;
                end
                           
            end % end of samples loop
            
            % store temporal vec (t) in cell
            allsubs_choicevec{1,subI}{1,id}   = t;
            
            clear t s 
            
        end % end of trials loop
    end % end of blocks loop 
    
    % update indx var so that it carries on after each participant
    counter                                     = counter + (indx/subI); 
    
end % end of subjects loop

% add data in one matrix
phase2_blockdata = [subj' blockno' trialno' numsamples' thisitem' thisprice' thisrank' reward' balance'];

% % save matrix in csv format for r and python
% csvwrite('economic_phase2_blockdata.csv', phase2_blockdata)

clear subj trialno blockno thisitem thisprice numsamples thisrank indx

%% EXTRACT ITEMS/PRICES AND RANKS %%

% loop over subjects
for subI = 1:nsubs
    
    tmpsub                  = find(phase2_blockdata(:,1) == subI);
    sub_blockdata           = phase2_blockdata((tmpsub),:);
    
    substruct.samples       = sub_blockdata(:,4);
    substruct.item          = sub_blockdata(:,5);
    substruct.price         = sub_blockdata(:,6);
    substruct.rank          = sub_blockdata(:,7);
    
    % store sub struct in cell
    allsubs_data{1,subI}    = substruct;
    
    clear tmpsub sub_blockdata
    
end % end of subject loop

%% STORE THE NUMBER OF SAMPLES AS A SEPARATE ARRAY FOR EEG

% for each subject average the number of samples to get one data value per
% subject. This array will then be used as a covariate in the individual
% differences analysis with EEG

for sub = 1:nsubs 

    % extract sub dataset
    tmp_sub                 = find(phase2_blockdata(:,1) == sub);
    sub_samples             = phase2_blockdata((tmp_sub),4);
    allsubs_samples(sub,1)  = nanmean(sub_samples);
    
    % clear stuff
    clear tmp_sub sub_samples 

end % end of subjects loop

save('avsamples.mat', 'allsubs_samples')

%% DEAL WITH SEQUENCES %%

% During the experiment I have been saving the index/items in the sequences
% but not the prices themselves. Below I create a new cell with sequences, only
% I store the prices (based on their index). This will be needed for model
% fitting/running ideal observer.

% further, create sequences with the ratings for prices/items in every
% sequence (i.e., link ratings with corresponding prices in every squence) 
% this will also be required for running the model

% loop over subjects
for sub = 1:nsubs
    
    % extract this_sub sequences
    sub_seq                                     = allsubs_sequences{1,sub};
    
    % extract this_sub indexes and prices
    sub_prices                                  = allsubs_prices{1,sub};
    
    % extract this_sub averaged ratings
    sub_rate                                    = allsubs_ratings{1,sub};
    
    % loop over sequences
    for seq = 1:size(sub_seq,2)
        
        % extract this_sequence
        this_seq                                = sub_seq{1,seq}';
        tmp_vec                                 = zeros(1,length(this_seq));
        tmp_rate                                = zeros(1,length(this_seq));
        
        % loop over itemns in this_seq
        for i = 1:length(this_seq)
            
            tmp_item                            = this_seq(1,i);
            tmp_vec(1,i)                        = sub_prices(tmp_item,2); % link this_seq items to their corresponding prices
            tmp_rate(1,i)                       = sub_rate(tmp_item,1); % link this_seq items to their corresponding ratings
        
        end % end of items loop
        
        % store new sequences 
        allsubs_price_sequences{1,sub}{1,seq}   = tmp_vec;
        allsubs_rate_sequences{1,sub}{1,seq}    = tmp_rate;
        
        % clear vars
        clear i tmp_item tmp_vec this_seq tmp_rate
        
    end % end of sequence loop
end % end of subjects loop

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
