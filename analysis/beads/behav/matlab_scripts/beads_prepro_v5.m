% PRE-PROCESSING AND ANALYSIS SCRIPT FOR THE BEADS TASK VERSION 5

% Part of the Optimal Stopping Problems Project

% Last update: 17/09/2022

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
% 4. average average sub draws and acc for each condition
% 5. run ideal observer 
% 6. average model draws and acc (total and for each condition)
% 7. run model fiiting 


%% INIT LOAD DATA %%

% clear workspace 
clear all
clc

% GET PATHS & DEFINE VARIABLES
% The four next lines (paths) should be changed to your paths 
startpath       = '/Users/christinadelta/githubstuff/rhul_stuff/SocialDeMa/';
modelfitpath    = fullfile(startpath, 'analysis', 'beads', 'behav', 'model_fitting');
iobserverpath   = fullfile(startpath, 'analysis', 'beads', 'behav', 'brunos_io');
resultspath     = fullfile(startpath, 'experiments', 'results');
croppedpath     = fullfile(startpath, 'analysis', 'beads', 'behav', 'cropped');

task            = 'beads';
subpath         = fullfile(resultspath, task);
session         = 1;

subs            = dir(fullfile(resultspath, task, '*sub*'));
nsubs           = length(subs);
% nsubs           = 2;

totaltrials     = 52; 
blocktrials     = 13;
maxdraws        = 10;
blocks          = 4;
conditions      = 2;
temp            = 0;
respoptions     = 3; % b,g,s

% init required vars
avdraws                 = nan(nsubs,1);
avacc                   = nan(nsubs,1);
easy_avdraws            = nan(nsubs,1);
diff_avdraws            = nan(nsubs,1);
easy_avacc              = nan(nsubs,1);
diff_avacc              = nan(nsubs,1);

% only keep subnames
subname                 = {subs.name};

for subI = 1:nsubs
        
%     if subI == 9 | subI == 10
%         continue
%     end
    
    fprintf('loading beads block data\n')  
    subject = subs(subI).name;
    subdir  = fullfile(resultspath, task,subject);
    fprintf('\t reading data from subject %d\n',subI); 
    
    co              = 1; % this will be used for spliting sequences and choice vectors in conditions one and two 
    ct              = 1; % this will be used for spliting sequences and choice vectors in conditions one and two 
    
    for blockI = 1:blocks
        
        fprintf('\t\t loading block %d\n\n',blockI);
        subFile = fullfile(subdir, sprintf('subject_%02d_task_%s_block_%02d_ses_%02d_logs.mat',subI, task, blockI,session));
        load(subFile)
        
        
        for trial = 1:blocktrials
             
            indx                            = ((blockI -1)*blocktrials) + trial; 
            
            block(indx)                     = blockI;
            trialno(indx)                   = logs.blocktrials(trial).trialnumber;
            urntype(indx)                   = logs.blocktrials(trial).urntype;
            
            % 1st draw is not stored so, add 1
            if logs.blocktrials(trial).draws < maxdraws
                draws(indx)                     = logs.blocktrials(trial).draws + 1;
            elseif logs.blocktrials(trial).draws == maxdraws
                draws(indx)                     = logs.blocktrials(trial).draws;
            end
            response(indx)                  = logs.blocktrials(trial).response;
            accuracy(indx)                  = logs.blocktrials(trial).accuracy;
            condition(indx)                 = logs.blocktrials(trial).condition;
            subj(indx)                      = subI;
            generaltrl(indx)                = indx;
            
            sequence                        = logs.blocktrials(trial).sequence; % extract tish trial sequence of draws
            
             % maybe here add this trial's responses. Given that I was not saving within sequence resposes but we need that info for
            % the model, I will create a nx3 matrix of choices for each sequence, when n=number of draws and 3=choice options (b,g,d)
            t                               = nan(draws(indx), respoptions); % init empty matrix
            
            for d = 1:draws(indx)
                if d ~= draws(indx) % if this is not the last draw add 0's to b and g columns and 1 to draw column
                    t(d,1:2)                = 0; % index zero for b and g columns
                    t(d,3)                  = 1; % index one for draw column
                else
                    if urntype(indx) == 1 & accuracy(indx) == 1 % this is a blue urn and sub was correct
                        t(d,2:3)            = 0; % index zero for g and draw columns
                        t(d,1)              = 1; % index one for b column (sub ressed blue)
                    elseif urntype(indx) == 1 & accuracy(indx) == 0 % this is a blue urn and sub was incorrect or did not respond
                        t(d,1)              = 0; % index zero for b 
                        t(d,2)              = 1; % index one for g column
                        t(d,3)              = 0; % index zero for draw 
                    elseif urntype(indx) == 0 & accuracy(indx) == 1 % this is a green urn and sub was correct
                        t(d,1)              = 0; % index zero for b 
                        t(d,2)              = 1; % index one for g column
                        t(d,3)              = 0; % index zero for s 
                    elseif urntype(indx) == 0 & accuracy(indx) == 0 % this is a green urn and sub was incorrect or did not respond
                        t(d,2:3)            = 0; % index zero for g and draw columns
                        t(d,1)              = 1; % index one for b column            
                    end
                end % end of if statement
            end % end of draws loop
            
            % add thistrial sequence and choive vector t in correct cell
            % based on condition
            if condition(indx) == 1 
                allsequences{1,subI}{1,condition(indx)}{1,co}       = sequence;
                allchoicevectors{1,subI}{1,condition(indx)}{1,co}   = t;
                
                co                                                  = co+1; % update co
            else 
                allsequences{1,subI}{1,condition(indx)}{1,ct}       = sequence;
                allchoicevectors{1,subI}{1,condition(indx)}{1,ct}   = t;
                
                ct                                                  = ct+1; % update co
            end
            
            clear t sequence    

        end % end of trial loop  
    end % end of block loop
    
    % add data in one matrix
    all_data = [subj' block' trialno' urntype' draws' response' accuracy' condition' generaltrl'];
    
    allsub_alldata{1,subI} = all_data;
    
    %% SPLIT SUB DATA INTO CONDITIONS
    
    for cond = 1:conditions % loop over conditions
        
        tmp                             = find(all_data(:,8) == cond);
        cond_data{subI}{cond}           = all_data((tmp),:);
        clear tmp
        
    end % end of condition loop
    
    %% AVERAGE PARTICIPANT DRAWS & ACCURACY %%
    
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
        
    end
    
    %% RUN IDEAL OBSERVER %%
    
    % first add modelpath to the path
    addpath(genpath(iobserverpath));
    
    % define parameters of the ideal observer
    R.alpha             = 1;            % softmax stochasticity parameter (for fitting to human behaviour) - this is not needed here
    R.error             = -10;          % cost for being wrong
    R.diff              = -20;          % The difference between the rewards for being correct (in this case no reward 10) and the cost of being wrong (-10).
    R.correct           = 10;           % reward for being correct
    R.q                 = [0.8 0.6];    % proportion of the majority value in sequence (60:40 split in this case)
    R.sample            = -0.25;            % the cost to sample
    
    for cond = 1:conditions
        
        thiscond_data    = cond_data{1,subI}{1,cond};
        thiscond_seq     = allsequences{1,subI}{1,cond};
        
        % what is the probability of this cond? 
        if cond == 1
            thisq = R.q(1);
        else 
            thisq = R.q(2);
        end
        
        R.thisq = thisq;
        
        % run backward induction
        [r, Qsat] = backWardInduction(thiscond_seq, R);
        
        % store ideal observer output
        io_output(cond).r       = r;
        io_output(cond).Qsat    = Qsat;
        
        % loop over condition trials to compute choices, picktrials and acc
        for i = 1: totaltrials/2
            
            choiceTrial             = find(squeeze(Qsat(i,:,3)) - max(squeeze(Qsat(i,:,1:2))') < 0); % which options this trial vec have an urn > sample
            pickTrial(i)            = choiceTrial(1); % pick the first of the choices
            [ma ma_i]               = max(squeeze(Qsat(i,pickTrial(i),:))); % which of the two urn was chosen? [based on AQ values]
            choice(i)               = ma_i;  % assign chosen urn
            choice(find(choice==2)) = 0; % recode it so it can be summed
            
        end
        
        % for each subject model instance and condition, calculate acc and
        % draws
        allsub_ioacc(subI,cond)     = mean(choice);
        allsubs_iodraws(subI,cond)  = mean(pickTrial);   
    
    end % end of condition loop
    
    % save this_sub ideal observer output
    allsubs_io{1,subI}              = io_output;
    
    clear R r Qsat
    
    %% RUN MODEL FITTING %%
%     % first add model fitting to the path
%     addpath(genpath(modelfitpath));
%     
%     % define parameters of the ideal observer
%     R.alpha             = 0.13;            % softmax stochasticity parameter (for fitting to human behaviour) - this is not needed here
%     R.error             = -10;          % cost for being wrong
%     R.diff              = -20;          % The difference between the rewards for being correct (in this case no reward 10) and the cost of being wrong (-10).
%     R.correct           = 10;           % reward for being correct
%     R.q                 = [0.8 0.6];    % proportion of the majority value in sequence (60:40 split in this case)
%     R.sample            = -0.25;            % the cost to sample
%     
%     for cond = 1:conditions
%         
%         % extract subject choice data, sequences
%         thisub_choices       = allchoicevectors{1,subI}{1,cond};
%         thisub_seq           = allsequences{1,subI}{1,cond};
%         cond_matrix          = cond_data{1,subI}{1,cond};
%        
%         % extract urn types form data matrix
%         info.urntypes        = cond_matrix(:,4);
%         info.condtrials      = totaltrials/conditions;
%         info.numdraws        = cond_matrix(:,5);
%         % Cs                   = -0.25;   
%         
%         % what is the probability of this cond? 
%         if cond == 1
%             thisq = R.q(1);
%         else 
%             thisq = R.q(2);
%         end
%         
%         R.thisq = thisq;
%         
%         
%         
%         
%     end % end of conditions loop
    
    
end % end of subject loop


