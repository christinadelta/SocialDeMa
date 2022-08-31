% PRE-PROCESSING AND ANALYSIS SCRIPT FOR THE BEADS TASK VERSION 4

% Part of the Optimal Stopping Problems Project

% Last update: 30/08/2022

%% IMPORTANT NOTE %%

% I save two different types of mat files. In 4 mat files I save the block
% information. This kind of file contains the trials of the given block, the
% number of draws for each trial, the balance of this trial (reward/loss),
% response, accuracy, etc..

% The rest mat files contain sequence/trial information, such as: each draw
% of the given trial, trial start, bead onset, response time, etc...

% TOTAL MAT FILES for each subject: 56
% BLOCK MAT FILES: 4 logs
% SEQUENCE MAT FILES: 52 logs

% I first extract block data, remove nans and save the data in a csv file 
% Then I extract the sequence data, remove nans and save the data in a csv file 

%% PREPROCESSING - ANALYSIS STEPS %%

% 1. extract sub block data and store based on conditions [logs, sequences]
% 2. store sub data based on conditions
% 3. average sub draws and acc 
% 4. average average sub draws and acc for each condition
% 5. run ideal observer 
% 6. average model draws and acc (total and for each condition)
% 7. run model fiiting 
% 8. get difference in AQ values
% 9. run stats (random effects (agent type) anova)

%% INIT LOAD DATA %%

% GET PATHS & DEFINE VARIABLES
% The four next lines (paths) should be changed to your paths 
startpath       = '/Users/christinadelta/githubstuff/rhul_stuff/SocialDeMa/';
modelfitpath    = fullfile(startpath, 'analysis', 'beads', 'behav', 'model_fitting');
iobserverpath   = fullfile(startpath, 'analysis', 'beads', 'behav', 'ideal_observer');
analysispath    = fullfile(startpath, 'analysis', 'beads', 'behav', 'prepro_data');
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
blocks          = 4;
conditions      = 2;
temp            = 0;
respoptions     = 3; % b,g,s
counter         = 1; 

% init required vars
avdraws                 = nan(nsubs,1);
avacc                   = nan(nsubs,1);
easy_avdraws            = nan(nsubs,1);
diff_avdraws            = nan(nsubs,1);
easy_avacc              = nan(nsubs,1);
diff_avacc              = nan(nsubs,1);

% only keep subnames
subname                 = {subs.name};


%% 1. EXTRACT AND SAVE THE BLOCK DATA %%

for subI = 1:nsubs
    
    if subI == 9
        continue
    end
       
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
            draws(indx)                     = logs.blocktrials(trial).draws;
            response(indx)                  = logs.blocktrials(trial).response;
            accuracy(indx)                  = logs.blocktrials(trial).accuracy;
            condition(indx)                 = logs.blocktrials(trial).condition;
            subj(indx)                      = subI;
            
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
            
            counter = counter + 1;
           
            clear t sequence    

        end % end of trial loop
         
    end % end of block loop
    
    % add data in one matrix
    all_data = [subj' block' trialno' urntype' draws' response' accuracy' condition'];
    
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
    
    % define parameters of the model
    alpha       = 10;            % softmax stochasticity parameter (for fitting to human behaviour)
    Cw          = -10;          % cost for being wrong
    Cd          = -20;          % The difference between the rewards for being correct (in this case no reward 10) and the cost of being wrong (-10).
    Cc          = 10;           % reward for being correct
    prob        = [0.8 0.6];    % proportion of the majority value in sequence (60:40 split in this case)
    Cs          = -0.25;        % the cost to sample
    
    for cond = 1:conditions
        
        thiscond_data    = cond_data{1,subI}{1,cond};
        thiscond_seq     = allsequences{1,subI}{1,cond};
        
        % what is the probability of this cond? 
        if cond == 1
            thisq = prob(1);
        else 
            thisq = prob(2);
        end
        
        % run ideal observer 
        [ll, pickTrial, dQvec, ddec, aQvec choice]  = estimateLikelihoodf(alpha,Cw,Cc,thisq,Cs,thiscond_seq,1);
        
        % save ideal observer output 
        io_output(cond).pickTrials                  = pickTrial;
        io_output(cond).dQvec                       = dQvec;
        io_output(cond).aQvec                       = aQvec;
        io_output(cond).choices                     = choice;
        io_output(cond).ddec                        = ddec;
        
        %% COMPUTE MEAN ACCURACY, DRAWS & POINTS %%
        
        % compute model accuracy 
        allsub_modelacc(subI,cond)                  = mean(choice == 1);
        
        % compute number of draws
        allsubs_modeldraws(subI,cond)               = mean(pickTrial);
        
        % compute points 
        % allsubs_modelpoints(subI,cond)              = (sum(choice == 1) * Cc) + (sum(choice == 2) * Cw) - (sum(pickTrial) * Cs);
        
    end % end of condition loop
    
    % store all sub model in one cell 
    allsubs_io{1,subI}                              = io_output;
    
    %% RUN MODEL FITTING %%
    
    % first add modelpath to the path
    addpath(genpath(modelfitpath));
    
    % define parameters of the model
    alpha                   = 1;            % softmax stochasticity parameter (for fitting to human behaviour)
    Cw                      = -10;          % cost for being wrong     
    Cc                      = 10;           % reward for being correct 
    cost_diff               = -20;          % The difference between the rewards for being correct (in this case no reward 0) and the cost of being wrong (-1000).
    q                       = [0.8 0.6];    % proportion of the majority value in sequence (60/40 split in this case)
    Cs                      = -0.25;        % the cost to sample
    aqvec_switch            = 1;            % still not sure why exactly this is needed 
    cnt                     = 1;            % counter
    
    for cond = 1:conditions
       
        thisub_choices       = allchoicevectors{1,subI}{1,cond};
        thisub_seq           = allsequences{1,subI}{1,cond};
        cond_matrix          = cond_data{1,subI}{1,cond};
       
        % extract urn types form data matrix
        info.urntypes        = cond_matrix(:,4);
        info.condtrials      = totaltrials/conditions;
        info.numdraws        = cond_matrix(:,5);
       
        % is this cond 0.8 or 0.6 probability? 
        if cond == 1
            prob        = q(1);
        else
            prob        = q(2);
        end
        
        info.p          = prob;
        
        % aaand fit the model 
        [mparams, lla, all_ll, aQvec] = bayesbeads(thisub_seq, thisub_choices, info, alpha, Cw, Cc, cost_diff, Cs, cond, subI);
        
        % store all model output
        model_output(cond).params   = mparams;
        model_output(cond).lla      = lla;
        model_output(cond).all_ll   = all_ll;
        model_output(cond).aQvec    = aQvec;
        
        % extract AQ values for each condition
        cond_aQ                     = aQvec;
        
        % extract AQ values for each sequence 
        for trial = 1:length(cond_aQ)
            
            trial_aQ    = cond_aQ{1,trial}
            
            % loop over draws (AQs for that trial) 
            for draw = 1:size(trial_aQ,1)
                
                thisAQ                      = trial_aQ(draw,:);
                allsub_dAQs{1,subI}(cnt,1)  = (thisAQ(1) - thisAQ(2)) - thisAQ(3)
                
                cnt                         = cnt+1 % update counter
            end % end of draws (AQs) loop
        end % end of trial loop
    end % end of conditiion loop
    
    % save model subject fitting output and dAQ
    % allsubs_dAQ{1,subI}     = dAQ;
    allsubs_model{1,subI}   = model_output;
    
    %% PREPARE FOR REGRESSION ANALYSIS OF CROPPED EEG WITH dAQ %%
    
    % get the sum of the total draws (should be the same with the total number of
    % EEG epochs)
    totaldraws              = sum(all_data(:,5));
   
    % count of the total draws starts after the presentation of the first
    % draws. However, epochs start WITH the first bead presentation, so we
    % need to add 52 to the total draws:
    totaldraws              = totaldraws + totaltrials;
    
    % load the croped eeg data
    eeg_tmp                 = load(fullfile(croppedpath, sprintf('cropped_data_sub_%02d.mat',subI)));
    eeg_data                = eeg_tmp.data;
    
    % split the eeg data in frontal and parietal 
    frontal_eeg             = eeg_data(1:9,:,:);
    parietal_eeg            = eeg_data(10:20,:,:);
    
    % average eeg samples and channels over trials (that would results in
    % one averaged data point for each trial/epoch)
    for trl = 1:totaldraws
        
        tmp_frontal                     = frontal_eeg(:,:,trl);
        tmp_parietal                    = parietal_eeg(:,:,trl);
        frontal_aveeg{1,subI}(trl,1)    = mean(tmp_frontal(:));  % average channels*samples for this trial
        parietal_aveeg{1,subI}(trl,1)   = mean(tmp_parietal(:)); % average channels*samples for this trial
        
        clear tmp_frontal tmp_parietal
    end % end of trl 
    

    clear all_data accuracy draws response urntype condition subj
    
end % end of subject loop 

%% JUST SAVE STUFF %%
% store average sub and model draws and acc in one matrix
allagent_avdraws(:,1) = easy_avdraws;
allagent_avdraws(:,2) = diff_avdraws;
allagent_avdraws(:,3) = allsubs_modeldraws(:,1);
allagent_avdraws(:,4) = allsubs_modeldraws(:,2);

allagent_avacc(:,1) = easy_avacc;
allagent_avacc(:,2) = diff_avacc;
allagent_avacc(:,3) = allsub_modelacc(:,1);
allagent_avacc(:,4) = allsub_modelacc(:,2);

% save needed matrices
save avdraws
save allsubs_io
save allsub_modelacc 
save allsubs_modeldraws 
save allsubs_modelpoints
save allsubs_dAQ
save allsubs_model
save frontal_aveeg
save parietal_aveeg
save allsub_dAQs

%% RUN WITHIN-SUBJECT REGRESSION ANALYSIS %%

% loop over subjects
for sub = 1:nsubs 
    
    if sub == 9 | sub == 10
        continue
    end
    
    sub_dAQ         = allsub_dAQs{1,sub};
    sub_pareeg      = parietal_aveeg{1,sub};
    sub_fronteeg    = frontal_aveeg{1,sub};
    
    % first I need to add nan values to the dAQ vector as we cannot run
    % regression with imbalanced data
    s               = length(sub_dAQ) + 1; % where to continue from?
    e               = length(sub_pareeg);
    sub_dAQ(s:e,1)  = nan;
    
    % run regression analysis using parietal_eeg data as the dependent and
    % dAQs as factor
    mdl = fitlm(sub_dAQ,sub_pareeg)
    parietal_betas(sub,1) = table2array(mdl.Coefficients(2,1));
    
    clear mdl
    
    mdl = fitlm(sub_dAQ,sub_fronteeg)
    frontal_betas(sub,1) = table2array(mdl.Coefficients(2,1));
    
end

%% RUN 1-SAMPLE T-TESTS %%

% run t-test on the parietal data 
[h,p,ci,stats] = ttest(parietal_betas)

[h,p,ci,stats] = ttest(frontal_betas)


