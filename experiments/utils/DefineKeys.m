function[set] = DefineKeys(taskNb, set)

% This subfunction is part of the "OPTIMAL STOPPING EXPERIMENTS". 
% It runs via the main script of each task using the task number (taskNb)

% it creates a list of keys for each experiment alongside with the "global" 
% keys used in all experiments

% GLOBAL KEYS
% 1. Escape key (allows subject to quit the experiment)
% 2. Space key (allows subject to start the trials after instructions)
KbName('UnifyKeyNames');
set.code20          = KbName('space');
set.code21          = KbName('ESCAPE');


if taskNb == 1 % if this is the  beads task
    %% TASK 1 KEYS
    KbName('UnifyKeyNames');
    set.code1       = KbName('1!'); % 1! = Blue
    set.code2       = KbName('2@'); % 2@ = Green
    set.code3       = KbName('3#'); % 3# = Draw again

    set.code7       = KbName('a'); % answer a
    set.code8       = KbName('b'); % answer b
    set.code9       = KbName('c'); % answer c
    set.code10      = KbName('d'); % answer d

elseif taskNb == 2 % if this is the the economic best-choice task
    % TASK 2 KEYS
    
    % Unpack variables
    phase               = set.phase;
    
    if phase == 2
        
        KbName('UnifyKeyNames');
        set.code1      = KbName('1!'); % choose option
        set.code2      = KbName('2@'); % continue sampling

    end % end of phase statement
    
elseif taskNb == 3 % if this is the facial attractiveness task
    %% TASK 3 KEYS
     
    % Unpack variables
    phase               = set.phase;
    
    if phase == 2
        
        KbName('UnifyKeyNames');
        set.code1       = KbName('1!'); % I would never choose that contract
        set.code2       = KbName('2@'); % 2@ = Green
         
    end % end of phase statement
    
end % end of if statement 

end