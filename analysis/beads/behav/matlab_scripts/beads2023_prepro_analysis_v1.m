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
%   - model 1: free parameter 1 = Cost-sample
%   - model 2: free parameter 1 = Cost-error
%   - model 3: free parameter 1 = beta
%   - model 4: free parameter 2 = Cost-sample, Cost-error
%   - model 5: free parameter 2 = Cost-sample, beta

% POTENTIAL ADDITIONAL MODELS TO INCLUDE (will decide after making the above model fit and parameter recovery work):
%   - model 6: free parameter 1 = discounting factor (gamma)
%   - model 7: free parameter 3 = Cost-sample, beta, discounting factor (gamma)

% 5. parameter recovery
% 6. model comparison
% 7. compute AQ differences [draw-again higher-option]
% 8. next analysis steps will be introduced soon...

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
bmodelfitpath   = fullfile(startpath, 'analysis', 'beads', 'behav', 'model_fit');
biobserverpath  = fullfile(startpath, 'analysis', 'beads', 'behav', 'ideal_observer');
resultspath     = fullfile(datapath, 'beads', 'behav');
croppedpath     = fullfile(startpath, 'analysis', 'beads', 'behav', 'cropped');
behavpath       = fullfile(startpath, 'analysis', 'beads', 'behav');
simpath         = fullfile(behavpath, 'simulated');
addpath(genpath(fullfile(behavpath, 'matlab_scripts'))); % add matlab_scripts to teh path
addpath(genpath(fullfile(simpath))); 

task            = 'beads';
subs            = dir(fullfile(resultspath, '*sub*'));
nsubs           = length(subs);
% nsubs           = 5;

totaltrials     = 52; 
conditions      = 2;
condtrials      = totaltrials/conditions;
nmodels         = 2;

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
all_acc             = [easy_avacc diff_avacc];
all_draws           = [easy_avdraws diff_avdraws];

% make a struct with all the vectors needed for the analysis
anova_struct        = struct('all_draws', all_draws, 'all_acc', all_acc,...
    'all_ioacc', all_ioacc, 'all_iodraws', all_iodraws);

% run anovas and pairwise comparisons 
output_struct_one   = runBehavStats(nsubs, anova_struct); % output will be used for plotting 

clear anova_struct

%% FIT MODELS %%

% first add model fitting to the path
addpath(genpath(bmodelfitpath));

% define fixed parameters used in all models 
R.error             = -10;          % cost for being wrong
R.correct           = 10;           % reward for being correct
R.q                 = [0.8 0.6];    % proportion of the majority value in sequence (60:40 split in this case)
R.difference        = -20;
R.Cs                = -0.25;
R.beta              = 3;

% how many models, which models?
model_names         = {'CostSample' 'CostError' 'CsCerror' 'beta' 'betaCs'};
model_num           = length(model_names);

% loop over models 
for model = 1:model_num

    if model == 1

        R.initsample   = R.Cs;
        R.freeparams   = 1;

    elseif model == 2

        R.initdiff     = R.difference;
        R.freeparams   = 1;

    elseif model == 3

        R.initsample   = R.Cs;
        R.initdiff     = R.difference;
        R.freeparams   = 2;

    elseif model == 4

        R.initbeta     = 3;
        R.freeparams   = 1;

     elseif model == 5

        R.initsample   = R.Cs;
        R.initbeta     = 3;
        R.freeparams   = 2;

    end

    R.model             = model;
     
    % loop over subjects 
    for sub = 1:nsubs
         
        subdata             = cond_data{1,sub};         % all data matrix
        sub_choicesvec      = allsub_choiceVec{1,sub};  % subject choices (two 26 by 3 matricies)
        sub_sequence        = allsub_sequences{1,sub};  % sequence that this-subject was presented with
        
        % loop over conditions
        for cond = 1:conditions

            % extract condition data
            sub_cond        = subdata{1,cond};
            cond_choices    = sub_choicesvec{1,cond};
            cond_sequence   = sub_sequence{1,cond};

            % what is the probability of this cond? 
            if cond == 1

                thisq       = R.q(1);
            else 
                thisq       = R.q(2);
            end
        
            R.thisq         = thisq;

            % extract sequences from cell and store in matrix (26x10)
            for s = 1:size(cond_sequence,2)
                thiscond_seqmat(s,:)                        = cond_sequence{1,s};
            end
            
            % what is the urntype?
            urntype                                         = sub_cond(:,4);

            for u = 1:length(urntype)
                if urntype(u) == 0 % if green urn switch index coding
                    seq_ones                                = find(thiscond_seqmat(u,:) == 1);
                    seq_twos                                = find(thiscond_seqmat(u,:) == 2);
                    thiscond_seqmat(u,seq_ones)             = 2;
                    thiscond_seqmat(u,seq_twos)             = 1;
                end 
            end

            % recode 2s to 0s for backward induction 
            thiscond_seqmat(find(thiscond_seqmat==2))       = 0;

            modeloutput         = fitAllModel(R,thiscond_seqmat,cond_choices,urntype);
            mout(cond).modelout = modeloutput;

            % for each subject and condition, store 

        end % end of conditions loop
    end % end of subjects loop
end % end of models loop

%% RUN STATISTICS ON BEHAVIOUR, IO & MODELS %%


%% RECOVER MODEL PARAMETERS %%

%% COMPARE & CHOOSE WINNING MODEL %%

%% PREDICT EEG RESPONSES (EVOKED AND TF) USING AQ DIFFERENCE VALUES %%

%% PLOT STUFF %% 
