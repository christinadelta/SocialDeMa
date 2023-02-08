% PRE-PROCESSING SCRIPT FOR BEST-CHOICE FACE TASK

% VERSION 2 (created February 2023)

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

clear all
clc

% get paths and define vars
startpath               = '/Volumes/DeepSpaceStuff/optimal_stopping_data/data/';
% resultspath             = fullfile(startpath, 'experiments', 'results');
task                    = 'facial_attractiveness';
datatype                = 'behav';
subpath                 = fullfile(startpath, task, datatype);
phase                   = 2; % change this when extracting phase 2 data
session                 = 2; % change this when extracting phase 2 data

subs                    = dir(fullfile(subpath, '*sub*'));
nsubs                   = length(subs);
% nsubs                   = 1;
phase1_blocks           = 20;
% phase1_blocktrials      = 40;
% phase1_totaltrials      = phase1_blocks*40;

phase2_blocktrials      = 20;
phase2_blocks           = 2;
respoptions             = 2; % accept vs decline
thisequence             = nan(1,phase2_blocktrials); % every sequence has different number of draws
temp                    = 0;

% only keep subnames (not sure if this will be used)
subname                 = {subs.name};

%% EXTRACT AND SAVE PHASE 2 BLOCK DATA %%

counter                 = 0; 
taskname                = 'face';

% loop over subjects
for sub = 1:nsubs

    fprintf('loading facial_attractiveness phase 2 block data\n')  
    subject = subs(sub).name;
    subdir  = fullfile(subpath,subject);
    fprintf('\t reading data from subject %d\n',sub); 

    % extract data 
    % loop over blocks
    for blockI = 1:phase2_blocks

        fprintf('\t\t loading block %d\n\n',blockI);
        subFile = fullfile(subdir, sprintf('subject_%02d_task_%s_block_%02d_ses_%02d_phase_%02d_blocktrials_logs.mat',sub,taskname,blockI,session,phase));
        load(subFile)
    
        blocktrials         = length(logs.blocktrials); % how many sequences?
    
        % loop over trials/sequences
        for trial = 1:blocktrials
                
            indx                                = counter + ((blockI -1)*blocktrials) + trial; 
            id                                  = ((blockI -1)*blocktrials) + trial;
            
            subj(indx)                          = sub;
            blockno(indx)                       = logs.blocktrials(trial).block;
            trialno(indx)                       = logs.blocktrials(trial).trialnumber;
            numsamples(indx)                    = logs.blocktrials(trial).numsamples;
            thisitem(indx)                      = logs.blocktrials(trial).chosenitem;
            
            allsubs_sequences{1,sub}{1,id}      = logs.blocktrials(trial).sequence;

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
            allsubs_choicevec{1,sub}{1,id}   = t;
            
            clear t s   
        end % end of trials loop

    end % end of blocks loop

    % update indx var so that it carries on after each participant
    counter                                  = counter + (indx/sub); 

end % end of subjects loop

% add data in one matrix
phase2_blockdata = [subj' blockno' trialno' numsamples' thisitem'];

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

save('avsamples_v2.mat', 'allsubs_samples')



