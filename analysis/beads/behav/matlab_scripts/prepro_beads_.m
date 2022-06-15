% PRE-PROCESSING SCRIPT FOR THE BEADS TASK

% Part of the Optimal Stopping Problems Project

%% IMPORTANT NOTE %%

% I save two different types of mat files. In 4 mat files I save the block
% information. This kind of file contains the trials of the given block, the
% number of draws for each trial, the balance of this trial (reward/loss),
% response, accuracy, etc..

% The rest mat files contain sequence/trial information, such as: each draw
% of the given trial, trial start, bead onset, response time, etc...

% TOTAL MAT FILES for each subject: 56
% BLOCK MAT FILES: 4
% SEQUENCE MAT FILES: 52

% I first extract block data, remove nans and save the data in a csv file 
% Then I extract the sequence data, remove nans and save the data in a csv file 

%% INIT LOAD DATA %%

% GET PATHS & DEFINE VARIABLES
% The four next lines (paths) should be changed to your paths 
startpath       = '/Users/christinadelta/githubstuff/rhul_stuff/SocialDeMa/';
modelfitpath    = fullfile(startpath, 'analysis', 'beads', 'behav', 'model_fitting');
modelpath       = fullfile(startpath, 'analysis', 'beads', 'behav', 'ideal_observer');
analysispath    = fullfile(startpath, 'analysis', 'beads', 'behav', 'prepro_data');
resultspath     = fullfile(startpath, 'experiments', 'results');

task            = 'beads';
subpath         = fullfile(resultspath, task);
session         = 1;

subs            = dir(fullfile(resultspath, task, '*sub*'));
nsubs           = length(subs);

totaltrials     = 52; 
blocktrials     = 13;
blocks          = 4;
conditions      = 2;
thisequence     = nan(1,blocktrials); % every sequence has different number of draws
temp            = 0;
respoptions     = 3; % b,g,s

% only keep subnames
subname         = {subs.name};

%% EXTRACT AND SAVE THE BLOCK DATA %%

for subI = 1:nsubs 
    
    fprintf('loading beads block data\n')  
    subject = subs(subI).name;
    subdir  = fullfile(resultspath, task,subject);
    fprintf('\t reading data from subject %d\n',subI); 
    
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
            rate(indx)                      = logs.blocktrials(trial).thisrate;
            condition(indx)                 = logs.blocktrials(trial).condition;
            balance(indx)                   = logs.blocktrials(trial).balance;
            subj(indx)                      = subI;
            generaltrial(indx)              = indx;
            
            % extract sequences from log file 
            sequences{subI, blockI, trial}  = logs.blocktrials(trial).sequence;
            
            % for each block and trial, create matricies with 3 columns and
            % rows of draw-length. This will correspond to responses.
            % Columns: [blue, green, sample] -- this is because for some
            % reason I did not save within-sequence responses!! 
            t                               = nan(draws(indx), respoptions);
            
            for d = 1:draws(indx)
                if d ~= draws(indx)
                    
                    t(d,1:2)                = 0; % index zero for b and g columns
                    t(d,3)                  = 1; % index one for s column
                else
                    if urntype(indx) == 1 & accuracy(indx) == 1
                        
                        t(d,2:3)            = 0; % index zero for g and w columns
                        t(d,1)              = 1; % index one for b column
                        
                    elseif urntype(indx) == 1 & accuracy(indx) == 0
                        
                        t(d,1)              = 0; % index zero for b 
                        t(d,2)              = 1; % index one for g column
                        t(d,3)              = 0; % index zero for s 
                    elseif urntype(indx) == 0 & accuracy(indx) == 1
                        
                        t(d,1)              = 0; % index zero for b 
                        t(d,2)              = 1; % index one for g column
                        t(d,3)              = 0; % index zero for s 
                        
                    elseif urntype(indx) == 0 & accuracy(indx) == 0
                        
                        t(d,2:3)            = 0; % index zero for g and s columns
                        t(d,1)              = 1; % index one for b column                      
                    end
                end
            end % end of for loop
            % add responses to 
            response_batches{subI, blockI, trial} = t;
            clear t
            
        end % end of trial loop
    end % end of block loop
end % end of subject loop

% add data in one matrix
all_data = [subj' block' trialno' urntype' draws' response' accuracy' rate' condition' balance'];

clear accuracy balance block urntype trialno draws response accuracy rate condition balance indx subj

% remove nans if any
% block_data(any(isnan(block_data), 2), :)  = []; (let's not remove nan's yet)

% save matrix in csv format in case we want to run analyses in r and/or python
csvwrite('beads_alldata.csv', all_data)

%% EXTRACT AND SAVE THE SEQUENCE DATA %%

% in the lines of code above, I extract info about block trials but the the
% information for each sequence (e.g. reaction time for deciding to either draw or choose an urn, each bead of the sequence...)  

% load subject
for subI = 1:nsubs
    
    fprintf('loading beads sequence data\n')  
    subject     = subs(subI).name;
    subdir      = fullfile(resultspath, task,subject);
    fprintf('\t reading data from subject %d\n',subI); 
    
    % load block
    for blockI = 1:blocks
        
        fprintf('\t\t loading block %d\n\n',blockI);
        
        % load trial
        for trial = 1:blocktrials
            
            % load sequence file
            fprintf('\t\t loading sequence %d\n\n',trial);
            subFile = fullfile(subdir, sprintf('subject_%02d_task_%s_block_%02d_trial_%02d_ses_%02d_draw_logs.mat',subI, task, blockI, trial,session));
            load(subFile)                     
            
            % how many draws in that sequence?
            thisequence(trial)          = length(logs.draws);
            
            for i = 1:thisequence(trial)
                
                index                   = temp + i;    

                thisblock(index)        = logs.draws(i).block;
                trialnb(index)          = logs.draws(i).trialnumber;
                drawno(index)           = logs.draws(i).thisdraw;
                bead(index)             = logs.draws(i).thisbead;
                rt(index)               = logs.draws(i).rt;
                subj(index)             = subI;
                
            end % end of sequence loop
            temp                        = temp + thisequence(trial); % update temp 
            
        end % end of trial loop
    end % end of block loop
end % end of subject loop

% add all sequence data in one matrix
sequence_data = [subj' thisblock' trialnb' drawno' bead' rt'];

% remove nans if any
% sequence_data(any(isnan(sequence_data), 2), :)  = []; (let's not remove nan's yet)

% save matrix in csv format in case we want to run analyses in r and/or python
csvwrite('beads_sequencedata.csv', sequence_data)

%% RUN MODEL FITTING %%

% first add modelpath to the path
addpath(genpath(modelfitpath));

% define parameters of the model
alpha           = 1;            % softmax stochasticity parameter (for fitting to human behaviour)
Cw              = -10;          % cost for being wrong     
cost_diff       = -20;          % The difference between the rewards for being correct (in this case no reward 0) and the cost of being wrong (-1000).
q               = [0.8 0.6];    % proportion of the majority value in sequence (60/40 split in this case)
Cs              = -0.25;        % the cost to sample
aqvec_switch    = 1;            % still not sure why exactly this is needed 

% loop over subjects 
for sub = 1:nsubs
    
    % loop over blocks 
    for block = 1:blocks
        
        % extract block data
        tmp             = find(all_data(:,2) == block);
        this_blockdata  = all_data((tmp),:);
        
        for trl = 1:blocktrials
            
            % init trial info struct
            trialinfo       = [];
            
            index           = ((block - 1) * blocktrials) + trl;
            
            % extract sub/block sequences
            this_sequence   = sequences{sub, block, trl};
            this_response   = response_batches{sub, block, trl};

            thisurn         = this_blockdata(trl,4); % blue or green urn?
            accurate        = this_blockdata(trl,7); % correct or incorrect?
            thiscond        = this_blockdata(trl,9); % 0.8 or 0.6 probabiltiy?
            
            % determine the probability (q) for this trial 
            if thiscond == 1
                thisq = q(1);
            else
                thisq = q(2);
            end
            
            % add trial info in one struct for modelling
            trialinfo.q         = thisq;
            trialinfo.urn       = thisurn;
            trialinfo.accurate  = accurate;
            trialinfo.cond      = thiscond;
            
            % model fitting 
            [mparams, lla, aQvec] = bayesbeads(this_sequence, this_response, trialinfo, alpha, Cw, cost_diff, Cs, sub);
            
            tempparams{sub}{thiscond}{index}   = mparams;
            templla{sub}{thiscond}{index}      = lla;
            tempaQvec{sub}{thiscond}{index}    = aQvec;

        end % end of trials loop
    end % end of block loop
end % end of sub loop


clear thisurn thisq thiscond accurate mparams lla aQvec index this_sequence this_response trialinfo sub

% remove empty cells from the model outputs 
for sub = 1:nsubs 
    for cond = 1:conditions
        modelparams{1,sub}{1,cond}  = tempparams{1,sub}{1,cond}(~cellfun('isempty',tempparams{1,sub}{1,cond}));
        modellla{1,sub}{1,cond}     = templla{1,sub}{1,cond}(~cellfun('isempty',templla{1,sub}{1,cond}));
        modelaQvec{1,sub}{1,cond}   = tempaQvec{1,sub}{1,cond}(~cellfun('isempty',tempaQvec{1,sub}{1,cond}));
        
    end 
end

%% RUN IDEAL OBSERVER %%

% first add modelpath to the path
addpath(genpath(modelpath));

% TODO:
% INstead of looping over blocks, loop and run ideal observer over conditions  

% define parameters of the model
alpha       = 1;            % softmax stochasticity parameter (for fitting to human behaviour)
Cw          = -10;          % cost for being wrong
Cd          = -20;          % The difference between the rewards for being correct (in this case no reward 10) and the cost of being wrong (-10).
Cc          = 10;           % cost for being correct
q           = [0.8 0.6];    % proportion of the majority value in sequence (60:40 split in this case)
cs          = -0.25;        % the cost to sample

% loop over subjects 
for sub = 1:nsubs
    
    % loop over blocks 
    for block = 1:blocks
        
        % extract block data
        tmp                     = find(block_data(:,2) == block);
        this_blockdata          = block_data((tmp),:);
        
        for trl = 1:blocktrials
            % extract sub/block sequences
            this_sequence       = sequences{sub, block, trl};

            thisurn             = this_blockdata(trl,4); % blue or green urn?
            accurate            = this_blockdata(trl,7); % correct or incorrect?
            thiscond            = this_blockdata(trl,9); % 0.8 or 0.6 probabiltiy?
            
            % determine the probability (q) for this trial 
            if thiscond == 1
                thisq           = q(1);
            else
                thisq           = q(2);
            end
        
            % run ideal observer 
            [ll, pickTrial, dQvec, ddec, aQvec choice] = estimateLikelihoodf(alpha,Cw,thisq,Cs,this_sequence,1);
            
            pick_trials{sub, block, trl}    = pickTrial;
            blockdQvec{sub, block, trl}     = dQvec;
            bloxkaQvec{sub, block, trl}     = aQvec;
            blockchoices(sub, block, trl)   = choice;
            blockddec{sub, block, trl}      = ddec;
  
        end % end of trial loop
    end % end of block loop
end % end of sub loop





