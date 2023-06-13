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
%   - model 1: free parameter 1 = beta
%   - model 2: free parameter 2 = Cost-sample, beta

% POTENTIAL ADDITIONAL MODELS TO INCLUDE (will decide after making the above model fit and parameter recovery work):
%   - model 7: free parameter 1 = discounting factor (gamma)
%   - model 8: free parameter 3 = Cost-sample, beta, discounting factor (gamma)

% 5. parameter recovery
% NOTE TO SELF:
% In parameter recovery, we first need to simulate 52 datasets and
% then simulate responses 
% TO TEST:
% To simulate responses I use backward induction and then the choice probabilities computed with the
% softmax function. Previously, I tried simulating responses using only the Q values but this seems incorrect. Not sure if using the softmax-choices will work, but it looks like the
% reasonable thing to do as a response model! -- Also, this is consistent with RL models 

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
%     R.diff              = -20;        % The difference between the rewards for being correct (in this case no reward 10) and the cost of being wrong (-10).
    R.correct           = 10;           % reward for being correct
    R.q                 = [0.8 0.6];    % proportion of the majority value in sequence (60:40 split in this case)
    R.sample            = -0.25;        % the cost to sample
    
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
model_names         = {'Beta' 'BetaCs'};
model_num           = length(model_names);

% loop over models 
for model = 1:model_num

    % define free parameters
    R.initsample    = R.Cs;
    R.initbeta      = R.beta;

    % how many free parameters?
    if model == 1
        R.freeparams    = 1;

    elseif model == 2

        R.freeparams    = 2;
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

            modeloutput                         = fitAllModel(R,thiscond_seqmat,cond_choices,urntype);
            % mout(cond).modelout               = modeloutput; % not needed for now

            % for each subject and condition, store most important outputs:
            allsubNLL(sub,cond)                 = modeloutput.NLL;
            
            if R.freeparams == 1
                allsubFitParams(sub,cond)       = modeloutput.fittedX; % beta 
                
            elseif R.freeparams == 2
                allsubFitParams.X1(sub,cond)    = modeloutput.fittedX(1); % Cs 
                allsubFitParams.X2(sub,cond)    = modeloutput.fittedX(2); % beta
               
            end

            allsubAvSamples(sub,cond)           = modeloutput.avSamples;
            allsubAvPerformance(sub,cond)       = modeloutput.modelPerformance;

        end % end of conditions loop
    end % end of subjects loop

    % store the above for each model 
    allModelsNLL{1,model}                   = allsubNLL;
    allModelsFitParams{1,model}             = allsubFitParams;
    allModelsAvSamples{1,model}             = allsubAvSamples;
    allModelsAvPerformance{1,model}         = allsubAvPerformance;

    clear allsubNLL allsubFitParams allsubAvSamples allsubAvPerformance
end % end of models loop

%% RUN STATISTICS ON BEHAVIOUR, IO & MODELS %%

% run stats for all model, human behaviour % io
% add draws and performance in one vec
all_acc                     = [easy_avacc diff_avacc];
all_draws                   = [easy_avdraws diff_avdraws];

% extract model samples and performance
% costSample_modelSamples     = allModelsAvSamples{1,1}; % cost-sample model
% CerrorReward_modelSamples   = allModelsAvSamples{1,2}; % cost-error & reward model
% CsCerrorReward_modelSamples = allModelsAvSamples{1,3}; % cost-sample, cost-error & reward model
% costDiff_modelSamples       = allModelsAvSamples{1,4}; % difference model
beta_modelSamples           = allModelsAvSamples{1,1}; % beta model
betaCs_modelSamples         = allModelsAvSamples{1,2}; % beta & cost-sample model

% costSample_modelPerf        = allModelsAvPerformance{1,1};
% CerrorReward_modelPerf      = allModelsAvPerformance{1,2};
% CsCerrorReward_modelPerf    = allModelsAvPerformance{1,3};
% costDiff_modelPerf          = allModelsAvPerformance{1,4};
beta_modelPerf              = allModelsAvPerformance{1,1};
betaCs_modelPerf            = allModelsAvPerformance{1,2};

% make a struct with all the vectors needed for the analysis
anova_struct        = struct('all_draws', all_draws, 'all_acc', all_acc,...
    'all_ioacc', all_ioacc, 'all_iodraws', all_iodraws,...
    'beta_modelSamples', beta_modelSamples, 'betaCs_modelSamples', betaCs_modelSamples,...
    'beta_modelPerf', beta_modelPerf, 'betaCs_modelPerf', betaCs_modelPerf);

% run stats function
output_struct_one   = runBehavStats(nsubs, anova_struct); % output will be used for plotting 

clear anova_struct

%% RECOVER MODEL PARAMETERS %%

% define some initial variables for the simulations
simvars.ntrials         = totaltrials;
simvars.maxDraws        = 10;
simvars.qvals           = [0.8 0.6];
simvars.conditions      = conditions;
simvars.contrials       = totaltrials / conditions;

% models to recover 
% simModels               = {'Csample' 'errorReward' 'CsErrorReward' 'beta' 'BetaCs'};
simModels               = {'Csample' 'beta' 'BetaCs'};
num_simModels           = length(simModels); 

for m = 1:num_simModels

    if m == 1
        
        % define parameters for simulations
        cs_bounds               = [-4 0]; % maximum and minimum cost to sample
        nbins                   = 16;
        % allCs                   = linspace(cs_bounds(1), cs_bounds(2), nbins+1);
        allCs                   = [-4 -3 -2 -1 -0.5 -0.25 0];
        simR.correct            = 10;
        simR.error              = -10;
        simR.difference         = -20;
        simR.initbeta           = 3;
        simR.freeparams         = 1;
        simR.model              = m;

        for thisCs = 1:length(allCs)

            simR.Cs             = allCs(thisCs); % use the
            simR.initsample     = allCs(thisCs);

            for cond = 1:conditions

                simR.cond               = cond;
                simoutput{1,cond}       = simBeadsData(simvars, simR);
                
                sim_drawsequence        = simoutput{1,cond}.simsequences;
                sim_choiceVecs          = simoutput{1,cond}.simchoicevec;
                sim_urntype             = simoutput{1,cond}.simurns;

                sim_modeloutput{1,cond}        = fitAllModel(simR,sim_drawsequence,sim_choiceVecs,sim_urntype);
                simX(thisCs,cond)              = simR.initsample;
                fitX(thisCs,cond)              = sim_modeloutput{1,cond}.fittedX;
                NLL(thisCs,cond)               = sim_modeloutput{1,cond}.NLL;
                fitSamples(thisCs,cond)        = sim_modeloutput{1,cond}.avSamples;
                fitPerf(thisCs,cond)           = sim_modeloutput{1,cond}.modelPerformance;

            end % end of conditions loop
        end % end of cs loop

    elseif m == 2 % if this is the beta model

        % define parameters for simulations
        beta_bounds             = [0 15]; % maximum and minimum cost to sample
        nbins                   = 16;
        % allbetas                   = linspace(beta_bounds(1), beta_bounds(2), nbins+1);
        allbetas                = [0 1 3 5 7 10 12 15];
        simR.correct            = 10;
        simR.error              = -10;
        simR.difference         = -20;
        simR.Cs                 = -0.25; % use the
        simR.initsample         = -0.25;
        simR.freeparams         = 1;
        simR.model              = m;

        for thisBeta = 1:length(allbetas)

            simR.beta = allbetas(thisBeta);
            simR.initbeta = allbetas(thisBeta);

            for cond = 1:conditions
                
                % simulate sequences and responses 
                simR.cond               = cond;
                simoutput{1,cond}       = simBeadsData(simvars, simR);
                
                % extract info
                sim_drawsequence        = simoutput{1,cond}.simsequences;
                sim_choiceVecs          = simoutput{1,cond}.simchoicevec;
                sim_urntype             = simoutput{1,cond}.simurns;
                
                % fit simulated data
                sim_modeloutput{1,cond}        = fitAllModel(simR,sim_drawsequence,sim_choiceVecs,sim_urntype);
                simX(thiBeta,cond)              = simR.initsample;
                fitX(thisBeta,cond)              = sim_modeloutput{1,cond}.fittedX;
                NLL(thisBeta,cond)               = sim_modeloutput{1,cond}.NLL;
                fitSamples(thisBeta,cond)        = sim_modeloutput{1,cond}.avSamples;
                fitPerf(thisBeta,cond)           = sim_modeloutput{1,cond}.modelPerformance;

            end % end of conditions loop
        end % end of betas loop

    elseif m == 3

        % define parameters for simulations
        cs_bounds               = [-4 0]; % maximum and minimum cost to sample
        beta_bounds             = [0 15]; % maximum and minimum cost to sample
        nbins                   = [16 16];
        allCs                   = [-4 -3 -2 -1 -0.5 -0.25 0];
        % allCs                   = linspace(cs_bounds(1), cs_bounds(2), nbins+1);
        % allbetas                   = linspace(beta_bounds(1), beta_bounds(2), nbins+1);
        allbetas                = [0 1 3 5 7 10 12 15];
        simR.correct            = 10;
        simR.error              = -10;
        simR.difference         = -20;
        simR.freeparams         = 2;
        simR.model              = m;

        for thisCs = 1:length(allCs)

            simR.Cs             = allCs(thisCs); % use the
            simR.initsample     = allCs(thisCs);

            for thisBeta = 1:length(allbetas)

                simR.beta       = allbetas(thisBeta);
                simR.initbeta   = allbetas(thisBeta);

                for cond = 1:conditions

                    % simulate sequences and responses 
                    simR.cond               = cond;
                    simoutput{1,cond}       = simBeadsData(simvars, simR);
                
                    % extract info
                    sim_drawsequence        = simoutput{1,cond}.simsequences;
                    sim_choiceVecs          = simoutput{1,cond}.simchoicevec;
                    sim_urntype             = simoutput{1,cond}.simurns;

                    % fit simulated data
                    sim_modeloutput{1,cond}             = fitAllModel(simR,sim_drawsequence,sim_choiceVecs,sim_urntype);

                    tmpsimX1(thisBeta,cond)              = simR.initsample;
                    tmpsimX2(thisBeta,cond)              = simR.initbeta;
                    tmpfitX1(thisBeta,cond)              = sim_modeloutput{1,cond}.fittedX(1);
                    tmpfitX2(thisBeta,cond)              = sim_modeloutput{1,cond}.fittedX(2);
                    tmpNLL(thisBeta,cond)                = sim_modeloutput{1,cond}.NLL;
                    tmpfitSamples(thisBeta,cond)         = sim_modeloutput{1,cond}.avSamples;
                    tmpfitPerf(thisBeta,cond)            = sim_modeloutput{1,cond}.modelPerformance;

                    tmpsimX.sample                      = tmpsimX1;
                    tmpsimX.beta                        = tmpsimX2;
                    tmpfitX.sample                      = tmpfitX1;
                    tmpfitX.beta                        = tmpfitX2;

                end % end of condition loop
            end % end of betas loop

            simX(thisCs).allsimX            = tmpsimX;
            fitX(thisCs).allfitX            = tmpfitX;
            NLL(thisCs).allNLL              = tmpNLL;
            fitSamples(thisCs).allSamples   = tmpfitSamples;
            fitPerf(thisCs).allPerf         = tmpfitPerf;
            
        end % end of Cs loop
    end % end of if statement

    % store all models recovered parameters 
    paramRec_simX{1,m}          = simX;
    paramRec_fitX{1,m}          = fitX;
    paramRec_NLL{1,m}           = NLL;
    paramRec_fitSamples{1,m}    = fitSamples;
    paramRec_fitPerf{1,m}       = fitPerf;
    
    clear simX fitX NLL fitSamples fitPerf sim_modeloutput simoutput
end % end of models loop

%% COMPARE & CHOOSE WINNING MODEL %%

%% PREDICT EEG RESPONSES (EVOKED AND TF) USING AQ DIFFERENCE VALUES %%

%% PLOT STUFF %% 
