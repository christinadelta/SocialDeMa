function[keys] = DefineKeys(taskNb)

% This subfunction is part of the "OPTIMAL STOPPING EXPERIMENTS". 
% It runs via the main script of each task using the task number (taskNb)

% it creates a list of keys for each experiment alongside with the "global" 
% keys used in all experiments:

taskNb = taskNb;

% GLOBAL KEYS
% 1. Escape key (allows subject to quit the experiment)
% 2. Space key (allows subject to start the trials after instructions)

KbName('UnifyKeyNames');
keys.code1      = KbName('1!'); % 1! = Blue
keys.code2      = KbName('2@'); % 2@ = Green
keys.code3      = KbName('3#'); % 3# = Draw again
keys.code4      = KbName('4$');
%keys.code4      = KbName('LeftArrow'); % not confident
% keys.code5      = KbName('DownArrow'); % moderately confident
% keys.code6      = KbName('RightArrow'); % very confident 
keys.code8      = KbName('space');
keys.code9      = KbName('ESCAPE');

% keys.code10     = KbName('a'); % answer a
% keys.code11     = KbName('b'); % answer b
% keys.code12     = KbName('c'); % answer c
% keys.code13     = KbName('d'); % answer d

end