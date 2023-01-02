% PRE-PROCESSING AND ANALYSIS SCRIPT FOR THE BEADS TASK VERSION 5

% Part of the Optimal Stopping Problems Project

% UPDATES:
% 30/09/2022
% 12/10/2022
% 18/10/2022

%%% Changes introduced: %%% 
% trying to run ideal observer using an adapted version of bruno's code
% trying to fit model using an adapted version of bruno's code
% This version only preprocesses behavioural data and runs the model.
% Regressions and correlations with the cropped EEG data are performed in
% beads_analysis_v2.m

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

%% PREPROCESSING STEPS %%

% 1. extract sub block data and store based on conditions [logs, sequences]
% 2. store sub data based on conditions
% 3. average sub draws and acc 
% 4. average sub draws and acc for each condition
% 5. run ideal observer 
% 6. average model draws and acc (total and for each condition)
% 7. run model fiiting 

% TODO:
% 1) Look at acc of both model-fitting outputs 


%% INIT LOAD DATA %%

% clear workspace 
clear all
clc

% GET PATHS & DEFINE VARIABLES
% The four next lines (paths) should be changed to your paths 
startpath       = '/Users/christinadelta/githubstuff/rhul_stuff/SocialDeMa/';
modelfitpath    = fullfile(startpath, 'analysis', 'beads', 'behav', 'model_fitting');
bmodelfitpath   = fullfile(startpath, 'analysis', 'beads', 'behav', 'brunos_modelfit');
iobserverpath   = fullfile(startpath, 'analysis', 'beads', 'behav', 'ideal_observer');
biobserverpath  = fullfile(startpath, 'analysis', 'beads', 'behav', 'brunos_io');
resultspath     = fullfile(startpath, 'experiments', 'results');
croppedpath     = fullfile(startpath, 'analysis', 'beads', 'behav', 'cropped');
behavpath       = fullfile(startpath, 'analysis', 'beads', 'behav');
addpath(genpath(fullfile(behavpath, 'matlab_scripts'))); % add matlab_scripts to teh path

task            = 'beads';
subs            = dir(fullfile(resultspath, task, '*sub*'));
nsubs           = length(subs);

totaltrials     = 52; 
conditions      = 2;

% init required vars
avdraws                 = nan(nsubs,1);
avacc                   = nan(nsubs,1);
easy_avdraws            = nan(nsubs,1);
diff_avdraws            = nan(nsubs,1);
easy_avacc              = nan(nsubs,1);
diff_avacc              = nan(nsubs,1);

for subI = 1:nsubs
    %% EXTRACT DATA FROM LOG FILES

    fprintf('loading beads block data\n')  
    subject = subs(subI).name;
    subdir  = fullfile(resultspath, task,subject);
    fprintf('\t reading data from subject %d\n',subI); 
    
    % extract blocktrial data
    [subsequences,subchoiceVec,all_data]    = get_blockdata(subdir,subI,task);
    
    % store in cell for each participant
    allsub_alldata{1,subI}                  = all_data;
    allsub_sequences{1,subI}                = subsequences;
    allsub_choiceVec{1,subI}                = subchoiceVec;
    
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
    
        % we will need: cost_correct, cost_wrong, cost_sample, numdraws, acc 
        allsub_points(subI, cond) = (sum(cond_acc==1)*10) + (sum(cond_acc==0)*-10) + (sum(cond_draws)*-0.25);
        
    end
    
    %% RUN IDEAL OBSERVER - BRUNO'S VERSION & NICK'S VERSION %%
    
    % first add modelpath to the path
    addpath(genpath(biobserverpath));
    addpath(genpath(iobserverpath));
    
    % define parameters of the ideal observer
    R.alpha             = 1;            % softmax stochasticity parameter (for fitting to human behaviour) - this is not needed here
    R.error             = -10;          % cost for being wrong
%     R.diff              = -20;          % The difference between the rewards for being correct (in this case no reward 10) and the cost of being wrong (-10).
    R.correct           = 10;           % reward for being correct
    R.q                 = [0.8 0.6];    % proportion of the majority value in sequence (60:40 split in this case)
    R.sample            = -0.25;            % the cost to sample
    R.samplen           = -0.025;
    
    for cond = 1:conditions
        
        thiscond_data    = cond_data{1,subI}{1,cond};
        thiscond_seq     = subsequences{1,cond};
        
        % extract sequences from cell and store in matrix (26x10)
        for i = 1:size(thiscond_seq,2)
            thiscond_seqmat(i,:) = thiscond_seq{1,i};
        
        end
        
        tmp_seqmat = thiscond_seqmat;
        
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
        thiscond_seqmat(find(thiscond_seqmat==2))=0;
        
        % run backward induction (bruno's code)
        [r, Qsat] = backWardInduction(thiscond_seqmat, R);
        
        % store ideal observer output
        bio_output(cond).r       = r;
        bio_output(cond).Qsat    = Qsat;
        
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
        allsub_bioacc(subI,cond)     = mean(choice);
        allsubs_biodraws(subI,cond)  = mean(pickTrial); 
        allsubs_biopoints(subI,cond) = (sum(choice==1)*R.correct) + (sum(choice==0)*R.error) + (sum(pickTrial)*R.sample);
        
        % run estimateLikelihoodf (Nick's code)
        [ll, picktrl, dQvec, ddec, aQvec choices] = estimateLikelihoodf_io(tmp_seqmat,R);
        
        io_output(cond).aQvec       = aQvec;
        io_output(cond).picktrl     = picktrl;
        
        choices(find(choices==2))   = 0; % re-code incorrect responses
        
        % store all io outpus
        allsub_ioacc(subI,cond)     = mean(choices==1);
        allsubs_iodraws(subI,cond)  = mean(picktrl); 
        allsubs_iopoints(subI,cond) = (sum(choices==1)*R.correct) + (sum(choices==0)*R.error) + (sum(picktrl)*R.sample);
 
        clear r Qsat thiscond_data thiscond_seq thiscond_seqmat choice pickTrial urntype
        
    end % end of condition loop
    
    % clear workspace 
    clear R cond cond_acc cond_draws urntype 
    
    % save this_sub ideal observer output
    allsubs_bio{1,subI}                 = bio_output; % bruno's io output
    allsubs_io{1,subI}                  = io_output; % Nick's io output
    
    
    %% RUN MODEL FITTING %%
    
    % for now (and for diagnostic purposes), will fit data using a) Nick's
    % model_fitting code and b) Bruno's code (will adapt it).
    % later on, we'll see which one will be used. I use Nick's model
    % fitting method to adapt it to Bruno's backwardInduction code
    
    % first add model fitting to the path
    addpath(genpath(bmodelfitpath));
    addpath(genpath(modelfitpath));
    
    % define parameters of the ideal observer
    R.initbeta          = 1;            % softmax stochasticity parameter (for fitting to human behaviour) - this is not needed here
    R.error             = -10;          % cost for being wrong
    R.diff              = -20;          % The difference between the rewards for being correct (in this case no reward 10) and the cost of being wrong (-10).
    R.correct           = 10;           % reward for being correct
    R.q                 = [0.8 0.6];    % proportion of the majority value in sequence (60:40 split in this case)
    R.initsample        = -0.25;        % the cost to sample
    R.initsamplen       = -0.025; 
    
    for cond = 1:conditions
        
        thiscond_data       = cond_data{1,subI}{1,cond};
        thiscond_seq        = subsequences{1,cond};
        thiscond_choiceVec  = subchoiceVec{1,cond};
        
        % what is the probability of this cond? 
        if cond == 1
            thisq = R.q(1);
        else 
            thisq = R.q(2);
        end
        
        R.thisq = thisq;
        
        % extract sequences from cell and store in matrix (26x10)
        for i = 1:size(thiscond_seq,2)
            thiscond_seqmat(i,:) = thiscond_seq{1,i};
        end
        
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
        thiscond_seqmat(find(thiscond_seqmat==2))=0;
        
        % fit free parameter using Nick's version of the model
        [mparams, lla, aQvec]           = bayesbeads(thiscond_seqmat, thiscond_choiceVec, R);
        
        % store model-fitting output
        allsubs_cs_mn(subI,cond)        = mparams(1);
        allsubs_beta_mn(subI,cond)      = mparams(2);
        allsubs_lla_mn(subI,cond)       = lla;
        allsubs_AQs_mn{1,subI}          = aQvec;
        
        % fit free parameter using Bruno's version of the model
        [minParams, ll, Qsad, cprob]    = bayesbeads_b(thiscond_seqmat, thiscond_choiceVec, R);
        
        allsubs_cs_mb(subI,cond)        = minParams(1);
        allsubs_beta_mb(subI,cond)      = minParams(2);
        allsubs_ll_mb(subI,cond)        = ll;
        allsubs_AQs_mb{1,subI}          = Qsad;
        
% %         % calculate model (fit) draws 
% %         for i = 1: totaltrials/2
% %             thist = squeeze(Qsad(i,:,:));
% %             for k = 1:length(thist)
% %                 
% %                 if any(isnan(thist(k,:)), 'all')
% %                     break
% %                 end
% %                 thisdraw = k;
% %             end
% %             
% %             choicet(i) = thisdraw -1; % to bring the draws back to normal (we don't want to include the last draw-choice)
% %             
% %             % remove nans from the model AQ/choice vector
% %             thist = rmmissing(thist);
% %             
% %             AQs_mb{1,i}= thist;
% %             
% %         end  
% %         allsubs_draws_mb(subI,cond) = mean(choicet);
% %         
% %         % store modelfit AQs in cell for both conditions 
% %         % this will be used for regressions with EEG data
% %         cond_AQs_mb{1,cond}  = AQs_mb;
        
    end % end of condition loop
    
    % store modelfit AQs in cell for all subs
    allsubs_AQs_mb{1,subI} = cond_AQs_mb;
    
       
end % end of subject loop


%% PLOT MODEL FIT RESULTS %%

allsub_draws(:,1) = easy_avdraws;
allsub_draws(:,2) = diff_avdraws;

% make plots
makePlots(allsubs_beta_mb, allsubs_cs_mb, allsubs_ll_mb, allsub_draws)
% makePlots(allsubs_cs_mb, allsubs_ll_mb, allsub_draws)

%% RUN BEHAV STATISTICS %%

% 1. 2x2 MIXED ANOVA on draws
% first add all the required data in one matrix 
subvec              = repmat(1:nsubs,1,4)';                             % create a vector with 4 copies participant number 
agentvec            = repmat([ones(1,nsubs*2) ones(1,nsubs*2)*2],1,1)'; % create a vector with 2 copies of agent type (indexed as 1=human, 2=io)
probvec             = repmat([ones(1,nsubs) ones(1,nsubs)*2],1,2)';     % create a vector with 2 copies of probability type (indexed as 1=0.8, 2=0.6)

% % add them all in one matrix for anova
% anovamat_draws(:,1) = subvec;
% anovamat_draws(:,2) = humanvec;
% anovamat_draws(:,3) = probvec;

% create 1 vec with all draws (human, io) 
drawsmat(:,1)       = easy_avdraws;
drawsmat(:,2)       = diff_avdraws;
drawsmat(:,3:4)     = allsubs_biodraws;
drawsvec            = drawsmat(:);

% create 1 vec with all acc (human, io) 
accmat(:,1)         = easy_avacc;
accmat(:,2)         = diff_avacc;
accmat(:,3:4)       = allsub_bioacc;
accvec              = accmat(:);

% run mixed 2x2 anova on draws 
[pvals,~,stats] = anovan(drawsvec, {subvec agentvec probvec}, ... 
'model','interaction', 'random',1,'varnames',{'subvec' 'agentvec' 'probvec'})

% run mixed 2x2 anova on accuracy 
[pvals,~,stats] = anovan(accvec, {subvec agentvec probvec}, ... 
'model','interaction', 'random',1,'varnames',{'subvec' 'agentvec' 'probvec'})



