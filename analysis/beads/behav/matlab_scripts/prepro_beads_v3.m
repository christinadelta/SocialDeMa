% PRE-PROCESSING AND ANALYSIS SCRIPT FOR THE BEADS TASK VERSION 3

% Part of the Optimal Stopping Problems Project

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

% 1. extract block data and store based on conditions [logs, sequences]
% 2. extract sequence data and store based on conditions 
% 3. Run ideal observer 
% 4. run model fiiting 

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
co              = 1; % this will be used for spliting sequences and choice vectors in conditions one and two
ct              = 1; % this will be used for spliting sequences and choice vectors in conditions one and two
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
end % end of subject loop

% add data in one matrix
all_data = [subj' block' trialno' urntype' draws' response' accuracy' rate' condition' balance'];

clear accuracy balance block urntype trialno draws response accuracy rate condition balance indx subj subI co ct d trial

% remove nans if any
% block_data(any(isnan(block_data), 2), :)  = []; (let's not remove nan's yet)

% save matrix in csv format in case we want to run analyses in r and/or python
% csvwrite('beads_alldata.csv', all_data)

%% RUN MODEL FITTING %%

% first split all_data into conditions the result should be a 26x10 matrix for
% each condition
for sub = 1:nsubs % loop over subjects 
    
    % extract this subject data 
    temp = find(all_data(:,1) == sub);
    sub_data = all_data((temp),:);
    
    for cond = 1:conditions % loop over conditions
        
        tmp                         = find(sub_data(:,9) == cond);
        cond_data{sub}{cond}        = sub_data((tmp),:);
        cond_data{sub}{cond}(:,11)  = 1:totaltrials/2; % just add row index number 
        clear tmp
    end % end of cond loop
end % end of subject loop

% define model parameters 
alpha                   = 1;            % softmax stochasticity parameter (for fitting to human behaviour)
Cw                      = -10;          % cost for being wrong     
cost_diff               = -20;          % The difference between the rewards for being correct (in this case no reward 0) and the cost of being wrong (-1000).
q                       = [0.8 0.6];    % proportion of the majority value in sequence (60/40 split in this case)
Cs                      = -0.25;        % the cost to sample
aqvec_switch            = 1;            % still not sure why exactly this is needed 

for sub = 1:nsubs
    
    % extract sub data, choices/responses, sequences 
    thisub_data         = cond_data{1,sub};
    thisub_choices      = allchoicevectors{1,sub};
    thisub_seq          = allsequences{1,sub};
    
    for cond = 1:conditions 
        % extarct condition data. choices, sequences
        cond_data       = thisub_data{1,cond};
        cond_choices    = thisub_choices{1,cond};
        cond_sequence   = thisub_seq{1,cond};
        
        % extract urn types form data matrix
        info.urntypes   = cond_data(:,4);
        info.condtrials = totaltrials/conditions;
        
        % is this cond 0.8 or 0.6 probability? 
        if cond == 1
            prob        = q(1);
        else
            prob        = q(2);
        end
        
        info.p          = prob;
        
        
        
        
    end
 
end % end of subjects loop 



