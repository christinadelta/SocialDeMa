% % PRE-PROCESSING AND ANALYSIS SCRIPT FOR THE BEADS TASK VERSION 1

% Part of the Optimal Stopping Problems Project

% CREATED: 21/03/2023
% This version will take over the previous (beads_prepro_v6.m) one which now is
% archived. 

% changes introduced in this version: 
% 1) I run everything using this script. The steps are described below 
% 2) I fixed the beta thingy after discussing with Nick 

%% PREPROCESSING STEPS %%

% 1. extract sub block data [logs, sequences]
% 2. store sub data based on conditions
% 3. average sub draws and acc for each condition

%% ANALYSES STEPS %%

% 1. run ideal observer (IO)
% 2. extract IO's sumpling rate and performance and average 
% 3. run statistical analyses using uaing the runBehavStats.m file:
%       a) run 2x2 anovas with factors: [probabilities, agent]
%       b) run pairwise comparisons to look at differences in model-human, 0.8 & 0.6 conditions

% 4. fit models:
%   model 1: one free params [beta]
%   model 2: two free params [beta & Cs]
%   parameter recovery
%   model comparison

% 5. compute AQ differences [draw-again higher-option]
% 6. next analysis steps will be introduced soon...

%% PLOTING %%

% 1. plot human behaviour and IO
% 2. plot human, IO and model fit behviours 
% 3. plot model fitting/parameter recovery stuff
% 4. plot regression scatterplots? -- for later 


%% IMPORTANT NOTE %%

% I save two different types of mat files. In 4 mat files I save the block
% information. This kind of file contains the trials of the given block, the
% number of draws for each trial, the balance of this trial (reward/loss),
% response, accuracy, etc..

% The rest mat files contain sequence/trial information, such as: each draw
% of the given trial, trial start, bead onset, response time, etc...

% TOTAL MAT FILES for each subject: 56
% BLOCK MAT FILES: 4 logs - USED for preprocessing
% SEQUENCE MAT FILES: 52 logs - NOT USED

% END OF PREAMBLE 

%% INIT LOAD DATA %%


% clear workspace 
clear all
clc

% GET PATHS & DEFINE VARIABLES
% The four next lines (paths) should be changed to your paths 
startpath       = '/Users/christinadelta/githubstuff/rhul_stuff/SocialDeMa/';
datapath        = '/Volumes/DeepSpaceStuff/optimal_stopping_data/data/';
bmodelfitpath   = fullfile(startpath, 'analysis', 'beads', 'behav', 'brunos_modelfit');
biobserverpath  = fullfile(startpath, 'analysis', 'beads', 'behav', 'brunos_io');
resultspath     = fullfile(datapath, 'beads', 'behav');
croppedpath     = fullfile(startpath, 'analysis', 'beads', 'behav', 'cropped');
behavpath       = fullfile(startpath, 'analysis', 'beads', 'behav');
addpath(genpath(fullfile(behavpath, 'matlab_scripts'))); % add matlab_scripts to teh path

task            = 'beads';
subs            = dir(fullfile(resultspath, '*sub*'));
% nsubs           = length(subs);
nsubs           = 5;

totaltrials     = 52; 
conditions      = 2;

% init required variables
avdraws                 = nan(nsubs,1);
avacc                   = nan(nsubs,1);
easy_avdraws            = nan(nsubs,1);
diff_avdraws            = nan(nsubs,1);
easy_avacc              = nan(nsubs,1);
diff_avacc              = nan(nsubs,1);

%% EXTRACT DATA FROM LOG FILES

for subI = 1:nsubs

   

    fprintf('loading beads block data\n')  
    subject = subs(subI).name;
    subdir  = fullfile(resultspath,subject);
    fprintf('\t reading data from subject %d\n',subI); 
    
    % extract blocktrial data
    [subsequences,subchoiceVec,all_data, draws_index]    = get_blockdata(subdir,subI,task);
    
    % store in cell for each participant
    allsub_alldata{1,subI}                  = all_data;
    allsub_sequences{1,subI}                = subsequences;
    allsub_choiceVec{1,subI}                = subchoiceVec;
    allsub_drawinfo{1,subI}                 = draws_index; % col1: draw, col2: trial

     %% SPLIT SUB DATA INTO CONDITIONS
    
    for cond = 1:conditions % loop over conditions
        
        tmp                             = find(all_data(:,8) == cond);
        cond_data{subI}{cond}           = all_data((tmp),:);
        clear tmp
        
    end % end of condition loop

    %% AVERAGE PARTICIPANT DRAWS & ACCURACY & CALCULATE POINTS%%
    
    % create a nx1 vector (n=number of participants) with the averaged number
    % of draws for each participant.
    % This vector will be used as a covariate for the individual differences
    % analysis in SPM12.
    
    sub_draws           = all_data(:,5);
    sub_acc             = all_data(:,7);
    avdraws(subI,1)     = mean(sub_draws);
    avacc(subI,1)       = mean(sub_acc);
    
    clear sub_draws sub_acc
    
    % average draws and acc for each condition
    for cond = 1:conditions
        
        tmp_cond                = cond_data{1,subI}{1,cond};
        cond_draws              = tmp_cond(:,5);
        cond_acc                = tmp_cond(:,7);
        
        if cond == 1
            easy_avdraws(subI,1) = mean(cond_draws);
            easy_avacc(subI,1)   = mean(cond_acc);
        else
            diff_avdraws(subI,1) = mean(cond_draws);
            diff_avacc(subI,1)   = mean(cond_acc);
            
        end
        
        % we will now calculate participant points. we can use this to se how each participant performed and to compare
        % each participant with their corresponding model instant 
    
        % we will need: cost_correct, cost_wrong, cost_to_sample, numdraws, acc 
        allsub_points(subI, cond) = (sum(cond_acc==1)*10) + (sum(cond_acc==0)*-10) + (sum(cond_draws)*-0.25);
        
    end

    %% RUN IDEAL OBSERVER %%
    
    % first add modelpath to the path
    addpath(genpath(biobserverpath));
    
    % define parameters of the ideal observer
    R.alpha             = 1;            % softmax stochasticity parameter (for fitting to human behaviour) - this is not needed here
    R.error             = -10;          % cost for being wrong
%     R.diff              = -20;          % The difference between the rewards for being correct (in this case no reward 10) and the cost of being wrong (-10).
    R.correct           = 10;           % reward for being correct
    R.q                 = [0.8 0.6];    % proportion of the majority value in sequence (60:40 split in this case)
    R.sample            = -0.25;            % the cost to sample
    
    for cond = 1:conditions
        
        thiscond_data               = cond_data{1,subI}{1,cond};
        thiscond_seq                = subsequences{1,cond};
        
        % extract sequences from cell and store in matrix (26x10)
        for i = 1:size(thiscond_seq,2)
            thiscond_seqmat(i,:)    = thiscond_seq{1,i};
        
        end
        
        % what is the probability of this cond? 
        if cond == 1
            thisq = R.q(1);
        else 
            thisq = R.q(2);
        end
        
        R.thisq = thisq;
        
        % what is the urntype?
        urntype = thiscond_data(:,4);
        
        for u = 1:length(urntype)
            if urntype(u) == 0 % if green urn switch index coding
                seq_ones = find(thiscond_seqmat(u,:) == 1);
                seq_twos = find(thiscond_seqmat(u,:) == 2);
                thiscond_seqmat(u,seq_ones) = 2;
                thiscond_seqmat(u,seq_twos) = 1;
            end 
        end
        
        % recode 2s to 0s for backward induction 
        thiscond_seqmat(find(thiscond_seqmat==2)) = 0;
        
        % run backward induction (bruno's code)
        [r, Qsat]                   = backWardInduction(thiscond_seqmat, R);
        
        % store ideal observer output
        io_output(cond).r           = r;
        io_output(cond).Qsat        = Qsat;
        
        % loop over condition trials to compute choices, picktrials and acc
        for i = 1: totaltrials/2
            
            choiceTrial             = find(squeeze(Qsat(i,:,3)) - max(squeeze(Qsat(i,:,1:2))') < 0); % which options this trial vec have an urn > sample
            pickTrial(i)            = choiceTrial(1); % pick the first of the choices
            [ma ma_i]               = max(squeeze(Qsat(i,pickTrial(i),:))); % which of the two urn was chosen? [based on AQ values]
            
            if (ma_i == 1 & urntype(i) == 1) | (ma_i == 2 & urntype(i) == 0)
                choice(i) = 1;
            else
                choice(i) = 0;
            end
        end  
        
        % for each subject model instance and condition, calculate acc and
        % draws
        all_ioacc(subI,cond)        = mean(choice);
        all_iodraws(subI,cond)      = mean(pickTrial); 
        all_iopoints(subI,cond)     = (sum(choice==1)*R.correct) + (sum(choice==0)*R.error) + (sum(pickTrial)*R.sample);
        
        clear r Qsat thiscond_data thiscond_seq thiscond_seqmat choice pickTrial urntype
        
    end % end of conditions loop
    
    % clear workspace 
    clear R cond cond_acc cond_draws urntype 
    
    % save this_sub ideal observer output
    all_io{1,subI}                 = io_output; % store entire io output


end % end of subjects loop

%% RUN STATISTICS ON BEHAVIOUR & IO %%

% add draws and performance in one vec
all_acc         = [easy_avacc diff_avacc];
all_draws       = [easy_avdraws diff_avdraws];


% run anovas and pairwise comparisons 
[draws_stats,acc_stats, pc_results, pc_tables] = ...
    runBehavStats(nsubs,all_draws,all_acc,all_ioacc,all_iodraws);



