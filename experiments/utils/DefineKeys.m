function[set] = DefineKeys(taskNb, set)

% This subfunction is part of the "OPTIMAL STOPPING EXPERIMENTS". 
% It runs via the main script of each task using the task number (taskNb)

% it creates a list of keys for each experiment alongside with the "global" 
% keys used in all experiments:

% Unpack variables
phase               = set.phase;

% GLOBAL KEYS
% 1. Escape key (allows subject to quit the experiment)
% 2. Space key (allows subject to start the trials after instructions)
KbName('UnifyKeyNames');
set.code20          = KbName('space');
set.code21          = KbName('ESCAPE');

if taskNb == 1 % if this is the  beads task

    set.code1       = KbName('1!'); % 1! = Blue
    set.code2       = KbName('2@'); % 2@ = Green
    set.code3       = KbName('3#'); % 3# = Draw again
    set.code4       = KbName('LeftArrow'); % not confident
    set.code5       = KbName('DownArrow'); % moderately confident
    set.code6       = KbName('RightArrow'); % very confident 

    set.code7       = KbName('a'); % answer a
    set.code8       = KbName('b'); % answer b
    set.code9       = KbName('c'); % answer c
    set.code10      = KbName('d'); % answer d
    
elseif taskNb == 2 % if this is the the economic best-choice task
    
    if phase == 1
        
        KbName('UnifyKeyNames');
        set.code1       = KbName('1!'); % I would never choose that contract
        set.code2       = KbName('2@'); % 2@ = Green
        set.code3       = KbName('3#'); % 3# = Draw again
        set.code4       = KbName('4$'); % not confident
        set.code5       = KbName('5%'); % moderately confident
        set.code6       = KbName('6^'); % very confident 
        set.code7       = KbName('7&'); % very confident 
        set.code8       = KbName('8*'); % answer a
        set.code9       = KbName('9('); % I would definitely choose that contract
        
    else % if this is phase 2
        
        set.code1      = KbName('1!'); % choose option
        set.code2      = KbName('2@'); % continue sampling
        
    end
    
end % end of if statement 

end