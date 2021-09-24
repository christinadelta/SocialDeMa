function [averaged] = GetAverage(taskName, wd, sub)

% LOAD AND PRE-PROCESS PHASE-ONE FACIAL-ATTRACTIVENESS DATA

% sub-function used only for the facial attractiveness task

% FOR EACH FACE, AVERAGE THE RATINGS AND SAVE THEM FOR USE IN THE SECOND
% FACE OF THE FACIAL ATTRACTIVENESS TASK

% 300 faces
% 2 ratings per face

%%  GET PATHS 
previous_phase  = 1;
resultsfolder   = fullfile(wd, 'results',taskName, sprintf('sub-%02d', sub));
session         = 1;
trials          = 20;
blocks          = 30;

%% LOAD THE BLOCK MAT FILES 

for block = 1:blocks
    
    % load the subfile
    subFile = fullfile(resultsfolder, sprintf('subject_%02d_task_%s_block_%02d_ses_%02d_phase_%02d_logs.mat',sub, taskName, block, session, previous_phase));
    load(subFile)
    
    for trial = 1:trials
        
        trial_index             = ((block - 1)*trials) + trial;
        
        % find trials with more than 1 rt and keep the 2nd 
        item(trial_index)       = logs.trials(trial).thisitem;
        rate(trial_index)       = logs.trials(trial).response;
        
    end % end of trial loop
end % end of block loop

%% STORE THE DATA IN ONE MATRIX AND GET THE AVERAGES 

data                = [item' rate'];

% how many items?
items               = length(unique(data(:,1)));
averaged            = zeros(1,items);

for i = 1:items
    
    % current item?
    thisitem        = data(:,1) == i;
    
    % for "thisitem" find the average of the two rates
    averaged(i)     = mean(data(thisitem,2));
        
end 

return