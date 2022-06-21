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
% csvwrite('beads_alldata.csv', all_data)

%% RUN MODEL FITTING %%

% first split all_data into conditions the result should be a 26x10 matrix for
% each condition
for sub = 1:nsubs % loop over subjects 
    
    for cond = 1:conditions % loop over conditions
        
        tmp                         = find(all_data(:,9) == cond);
        cond_data{sub}{cond}        = all_data((tmp),:);
        cond_data{sub}{cond}(:,11)  = 1:totaltrials/2; % just add row index number 
        clear tmp
    end % end of cond loop
end % end of subject loop

% split sequences and behvavioural responses into conditions 
for subi = 1:nsubs
    % extract subject data
    tmp = find(all_data(:,1) == subi);
    subdata = all_data((tmp),:);
    
    for block = 1:blocks
        % extract block data
        tmpblock = find(subdata(:,2) == block);
        blockdata = subdata((tmpblock),:);
        
        for trial = 1:blocktrials
            
            index = ((block -1)*blocktrials) + trial;
            if blockdata(trial,9) == 1 % if easy condition
                
                condsequence{1,subi}{1,1}{1,index} = sequences{subi, block, trial};
                condresponse{1,subi}{1,1}{1,index} = response_batches{subi, block, trial};
            else
                condsequence{1,subi}{1,2}{1,index} = sequences{subi, block, trial};
                condresponse{1,subi}{1,2}{1,index} = response_batches{subi, block, trial};
                
            end
        end % end of trial loop
    end % end of blocks loop   
end % end of subject loop

% clear empty cell arrays
for sub = 1:nsubs 
    for cond = 1:conditions
        condseq{1,sub}{1,cond} = condsequence{1,sub}{1,cond}(~cellfun('isempty',condsequence{1,sub}{1,cond}));
        condresp{1,sub}{1,cond} = condresponse{1,sub}{1,cond}(~cellfun('isempty',condresponse{1,sub}{1,cond}));
        
    end
    
end

clear tmp subdata tmpblock blockdata condsequence condresponse sub subi cond block trial

% ok, now it's time to run model fitting!!
% first add modelpath to the path
addpath(genpath(modelfitpath));

% define parameters of the model
alpha           = 1;            % softmax stochasticity parameter (for fitting to human behaviour)
Cw              = -10;          % cost for being wrong     
cost_diff       = -20;          % The difference between the rewards for being correct (in this case no reward 0) and the cost of being wrong (-1000).
q               = [0.8 0.6];    % proportion of the majority value in sequence (60/40 split in this case)
Cs              = -0.25;        % the cost to sample
aqvec_switch    = 1;            % still not sure why exactly this is needed 

% loop over subjects and conditions and run the model for all sequences
% within each condition
for sub = 1:nsubs
    
    
    
    
    
end



