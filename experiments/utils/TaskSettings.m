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
if taskNb == 1 % if this is the beads task
   
    % BASIC IMPORTANT SETTINGS
    set.blocks      = 4;    % number of blocks 
    % set.blinks    = 3;    % run only with eeg - instruct subject to blink after n number of trials 

    % EXPERIMENTAL SETTINGS
    set.welcomedur  = 2.5;  % welcome screen duration = 2.5 sec
    set.infoscreen  = 2.5;  % this screen appears at the beginning of every sequence and informs the participant abou sequence number and probabilities of draws 
    set.bead_dur    = 0.6;  % bead duration in seconds
    set.response    = 2.5;  % self-paced or up to 2.5 sec 
    set.confrating  = 15;   % duration of the confidence rating screen. Self-paced or up to 15 sec
    set.fix_dur     = 0.5;  % duration of the fixation cross
    set.feed_dur    = 2.5;  % duration of the feedback window self-paced or up to 3 sec
    set.isi         = .6;   % in seconds
    set.jitter      = .5;   % 4 sec

    % TASK PARAMETERS 
    set.trials      = 52;                       % total trials
    set.blocktrials = set.trials/set.blocks;    % number of trials per block
    set.draws       = 10;                       % draws per sequence/trial
    set.conds       = 2;
    set.prob        = [0.8 0.6];                % 8:2 or 6:4 proportion of the beads in each of the two urns
    set.penalty     = 0.25;                     % every time subject chooses to draw they are panished with a £0.25 loss
    set.win         = 10;                       % £10 pounds if the subject gets the urn right
    set.balance     = 0;                        % balance starts from zero
    
end % end of if statement 

end 
    