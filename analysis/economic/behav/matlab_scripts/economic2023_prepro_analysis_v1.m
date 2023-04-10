% COMPLETE PREPROCESSING % ANALYSIS OF ECONOMIC BEHAVIOURAL DATA VERSION 1 
% CREATED: APRIL 2023 

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

clc 
clear all

%% INIT LOAD DATA %%

% get paths and define vars
datapath               = '/Volumes/DeepSpaceStuff/optimal_stopping_data/data/';
startpath               = '/Users/christinadelta/gitrepos/SocialDeMa/';
behavpath               = fullfile(startpath, 'analysis', 'economic', 'behav');
task                    = 'economic';
datatype                = 'behav'; 
subpath                 = fullfile(datapath, task, datatype);
phase                   = 1; % change this when extracting phase 2 data
session                 = 1; % change this when extracting phase 2 data

subs                    = dir(fullfile(subpath, '*sub*'));
% nsubs                   = length(subs);
nsubs                   = 7;

% not sure which ones are needed but we'll leave them for now
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

% add directories to the path
addpath(fullfile(behavpath, 'matlab_scripts')); % add matlab_scripts to teh path

%% EXTRACT AND SAVE PHASE 1 DATA %%

% loop over subs
for sub = 1:nsubs

    fprintf('loading economic phase 1 data\n')  
    subject             = subs(sub).name;
    subdir              = fullfile(subpath,subject);
    fprintf('\t reading data from subject %d\n',sub); 

    [phase1_data,meanRatings]           = get_meanRatings(subdir,sub,task);

    % store each subjects mean ratings 
    all_ratings{1,sub}                  = meanRatings;
    all_phase1data{1,sub}               = phase1_data;

end % end of subjects loop

%% EXTRACT AND SAVE PHASE 2 (SEQUENCE) DATA %%

% During the experiment I have been saving the index/items in the sequences
% but not the prices themselves. Below I create a new cell with sequences, only
% I store the prices (based on their index). This will be needed for model
% fitting/running ideal observer.

% further, create sequences with the ratings for prices/items in every
% sequence (i.e., link ratings with corresponding prices in every squence) 
% this will also be required for running the model

% loop over subs
for sub = 1:nsubs

    fprintf('loading economic phase 2 block data\n')  
    subject = subs(sub).name;
    subdir  = fullfile(subpath,subject);
    fprintf('\t reading data from subject %d\n',sub); 

    [sequences,choicevecs,phase2_data, samples] = get_blockdata(subdir,sub,task);

    allsub_sequences{1,sub}                     = sequences;
    allsub_choiceVecs{1,sub}                    = choicevecs;
    allsub_samples(:,sub)                       = samples;

    % deal with subject data
    allsub_meanssamples(sub,1)                  = nanmean(samples); % for eeg

end % end of subjects loop

%% DEAL WITH SEQUENCES %%

% loop over subjects
for sub = 1:nsubs

    
    sub_seq         = allsub_sequences{1,sub}; % extract this_sub sequences
    sub_prices      = all_ratings{1,sub}(:,3); % extract this_sub indexes and prices
    sub_rate        = all_ratings{1,sub}(:,2); % extract this_sub averaged ratings

   % loop over sequences
    for seq = 1:size(sub_seq,2)
        
        % extract this_sequence
        this_seq                                = sub_seq{1,seq}';
        tmp_vec                                 = zeros(1,length(this_seq));
        tmp_rate                                = zeros(1,length(this_seq));
        
        % loop over itemns in this_seq
        for i = 1:length(this_seq)
            
            tmp_item                            = this_seq(1,i);
            tmp_vec(1,i)                        = sub_prices(tmp_item,1); % link this_seq items to their corresponding prices
            tmp_rate(1,i)                       = sub_rate(tmp_item,1); % link this_seq items to their corresponding ratings
        
        end % end of items loop
        
        % store new sequences 
        allsub_price_sequences{1,sub}{1,seq}   = tmp_vec;
        allsub_rate_sequences{1,sub}{1,seq}    = tmp_rate;
        
        % clear vars
        clear i tmp_item tmp_vec this_seq tmp_rate
        
    end % end of sequence loop
end % end of sequence loop
