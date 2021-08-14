function [set] = TaskSettings(taskNb)

% THIS IS A SUBFUNCTION, PART OF THE "OPTIMAL STOPPING EXPERIMENTS". 

% it takes as input the task number (from the main script) and outputs a
% structure with a few important startup parameters for each
% task

% NOTE:
% All parameters in this function can change as appropreate 

% GLOBAL PARAMETERS:
set.taskNb      = taskNb;   % initialize settings structure
set.fixation    = '+';      % fixation cross 
set.EEG         = 0;        % set to 1 when running in the eeglab

% create a list of settings and parameters for the rts task 
% BASIC IMPORTANT SETTINGS
set.blocks      = 3;    % number of blocks 
% set.blinks    = 3;    % run only with eeg - instruct subject to blink after n number of trials 

% EXPERIMENTAL SETTINGS
set.instr_dur   = 2.5;  % instruction duration = 2.5 sec
set.bead_dur    = 1;    % bead duration in seconds
set.response    = 2.5;  % duration of theresponse window 
set.feed_dur    = 3;    % duration of the feedback window in sec
set.isi         = 1;    % in seconds
set.jitter      = 2;    % 4 sec

% TASK PARAMETERS 
set.trials      = 48;                       % total trials
set.blocktrials = set.trials/set.blocks;    % number of trials per block
set.draws       = 10;                       % draws per sequence/trial
set.conds       = 2;
set.prob        = [0.8 0.6];                % 8:2 or 6:4 proportion of the beads in each of the two urns

end 
    