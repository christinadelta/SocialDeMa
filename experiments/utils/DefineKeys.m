function[keys] = DefineKeys(taskNb)

% This subfunction is part of the "OPTIMAL STOPPING EXPERIMENTS". 
% It runs via the main script of each task using the task number (taskNb)

% it creates a list of keys for each experiment alongside with the "global" 
% keys used in all experiments:

keys.taskNb = taskNb;

% GLOBAL KEYS
% 1. Escape key (allows subject to quit the experiment)
% 2. Space key (allows subject to start the trials after instructions)

KbName('UnifyKeyNames');
keys.code1        = KbName('1!'); % z animate
keys.code2        = KbName('2@'); % m inanimate
keys.code3        = KbName('3#');
keys.code8        = KbName('space');
keys.code9        = KbName('ESCAPE');

end