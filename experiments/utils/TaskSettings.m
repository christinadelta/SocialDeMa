function [set] = TaskSettings(taskNb, sess)

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
set.welcomedur  = 2.5;      % welcome screen duration = 2.5 sec
set.jitter      = .4;       % 0.4 sec
set.isi         = .4;       % in seconds

 % create a list of settings and parameters for the rts task 
if taskNb == 1 % if this is the beads task
   
    % BASIC IMPORTANT SETTINGS
    set.blocks      = 4;    % number of blocks 

    % EXPERIMENTAL SETTINGS
    set.infoscreen  = 2.5;  % this screen appears at the beginning of every sequence and informs the participant abou sequence number and probabilities of draws 
    set.bead_dur    = 0.5;  % bead duration in seconds
    set.response    = 2.5;  % self-paced or up to 2.5 sec 
    set.fix_dur     = 0.5;  % duration of the fixation cross
    set.feed_dur    = 1;    % duration of the feedback window self-paced or up to 3 sec

    % TASK PARAMETERS 
    set.trials      = 52;                       % total trials
    set.blocktrials = set.trials/set.blocks;    % number of trials per block
    set.draws       = 10;                       % draws per sequence/trial
    set.conds       = 2;
    set.prob        = [0.8 0.6];                % 8:2 or 6:4 proportion of the beads in each of the two urns
    set.penalty     = 0.25;                     % every time subject chooses to draw they are panished with a £0.25 loss
    set.win         = 10;                       % £10 pounds if the subject gets the urn right
    set.loss        = 10;                       % £10 pounds if the subject gets the urn wrong
    set.balance     = 0;                        % balance starts from zero
    
    % DEFINE EEG TRIGGER CODES IF EEG = 1
    if set.EEG == 1
        
        % 1. start with the sequence related triggers
        set.trigger1    = 1;    % BLUE URN - blue bead (high prob)
        set.trigger2    = 2;    % BLUE URN - green bead (low prob)
        set.trigger3    = 3;    % GREEN URN - green bead (high prob)
        set.trigger4    = 4;    % GREEN URN - blue bead (low prob)
        
        set.trigger5    = 5;    % response prompt
        set.trigger6    = 6;    % blue urn (choice)
        set.trigger7    = 7;    % green urn (choice)
        set.trigger8    = 8;    % draw again
        
        % 2. Confidence rating and feedback relatd triggers
        set.trigger9    = 9;    % confidence screen
        set.trigger10   = 10;   % rating 1 (not confident)
        set.trigger11   = 11;   % rating 2 (a little confident)
        set.trigger12   = 12;   % rating 3 (confident)
        set.trigger13   = 13;   % rating 4 (very confident)
        
        set.trigger14   = 14;   % feedback screen (you win!)
        set.trigger15   = 15;   % feedback screen (you lose!)
        set.trigger16   = 16;   % feedback screen (you lose - out of draws)
        set.trigger19   = 17;   % feedback (didn't respond)
        
        % 3. main script triggers 
        set.trigger100  = 100;  % condition trigger -- condition (easy)
        set.trigger101  = 101;  % condition trigger -- condition (difficult) 
        set.trigger102  = 102;  % sequence start 
        set.trigger103  = 103;  % sequence end
    end % end of if EEG statement
    
elseif taskNb == 2
    
    phase           = sess;     % this is only for the economic and facial attarctivenes tasks
    set.phase       = phase;

    if phase == 1
        
        % EXPERIMENTAL SETTINGS
        set.fix_dur         = .7;   % in sec
        
        % TASK/PHASE SETTINGS
        set.blocks          = 120; 
        set.itemReps        = 2; 
        set.totaltrials     = 480*set.itemReps;
        set.blocktrials     = set.totaltrials/set.blocks;
        
    else % if phase is 2
        
        % EXPERIMENTAL SETTINGS
        set.fix_dur         = .5;   % in sec
        set.stimdur         = 1;  % in sec
        set.response        = 2.5;  % indicative of max response time 
        set.feedback        = 3;    % this will be just the presentation of the accepted contract
        set.reward          = [1 0.5 0.25]; % rewards best on the 3 ranks
        set.balance         = 0;
        
        % TASK/PHASE SETTINGS
        set.blocks          = 3; 
        set.samples         = 10; 
        set.totaltrials     = 30;
        set.blocktrials     = set.totaltrials/set.blocks;
        set.phaseitems      = set.samples*set.totaltrials; % that's ~62% of the total phase 1 contracts
        
        % DEFINE EEG TRIGGERS 
        if set.EEG == 1
            
            set.trigger11   = 11; % response trigger -- subject accepted a contract
            set.trigger12   = 12; % response trigger -- subject samples again
            set.trigger13   = 13; % feedback trigger
            
            set.trigger100  = 100; % sequence start
            set.trigger101  = 101; % sequence end
        end
       
    end % end of phase if statement
    
elseif taskNb == 3
    
    phase                   = sess;         % this is only for the economic and facial attarctivenes tasks
    set.phase               = phase;
    
    if phase == 1
        
        % EXPERIMENTAL SETTINGS
        set.fix_dur         = .7;   % in sec
        set.response        = 10;   % indicative of max response time 
        
        % STIMULUS SETTINGS
        set.stimsize        = 250;  % resize images or not?
        set.stimsize_deg    = 4;    % degrees of visual angle
        
        % TASK/PHASE SETTINGS
        set.blocks          = 120; 
        set.itemReps        = 2; 
        set.totaltrials     = 480*set.itemReps;
        set.blocktrials     = set.totaltrials/set.blocks;
        
    else % if phase is 2
        
        % EXPERIMENTAL SETTINGS
        set.fix_dur         = .5;   % in sec
        set.stimdur         = 1.5;    % in sec
        set.response        = 2.5;  % indicative of max response time 
        set.feedback        = 3;    % this will be just the presentation of the accepted contract
        
        
        % STIMULUS SETTINGS
        set.stimsize        = 200;  % resize images or not?
        set.stimsize_deg    = 4;    % degrees of visual angle
        
        % TASK/PHASE SETTINGS
        set.blocks          = 3; 
        set.samples         = 10; 
        set.totaltrials     = 30;
        set.blocktrials     = set.totaltrials/set.blocks;
        set.phaseitems      = set.samples*set.totaltrials; % that's ~62% of the total phase 1 faces
        
        % DEFINE EEG TRIGGERS 
        if set.EEG == 1
            
            set.trigger11   = 11; % response trigger -- subject accepted a contract
            set.trigger12   = 12; % response trigger -- subject samples again
            set.trigger13   = 13; % feedback trigger
            
            set.trigger100  = 100; % sequence start
            set.trigger101  = 101; % sequence end
        end

    end % end of phase statement 
    
end % end of if statement 

end 
    