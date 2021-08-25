%% STARTUP SCRIPT 

% CLEAN UP        
clear;
clc  
close all hidden;

% THE SCRIPT SHOULD RUN IN DIRECTORY : '/SocialDeMa'
% PRESS RUN OR TYPE "startup_script" IN THE COMMAND WINDOW AND ENTER 
% THE TASK CODE NAME. THIS WILL ADD THE CORRECT PATHS OF THE CURRENT
% EXPERIMENT TO YOUR MATLAB PATH, THEN IT WILL RUN THE CORRECT "main_task"
% SCRIPT.

% TASK CODE NAMES:

% beads             = beads task 
% economic          = economic best-choice task
% face              = facial attractiveness best-choice task

% WORKING DIRECTORIES AND INFO:

% startpath (or root)                   = /SocialDeMa/
% working_dir                           = /startpath/experiments

% core experimental functions           = /working_dir/utils
% experimental stimuli                  = /working_dir/stimuli
% tasks directory                       = /working_dir/tasks
% participant log files                 = /working_dir/results


% DEFINE INITIAL PATHS
startpath           = pwd;
workingpath         = fullfile(startpath, 'experiments', 'tasks');
% taskpath            = fullfile(workingpath, 'tasks');

% create a user input dialog to gather information
prompt          = {'Enter task name (e.g. beads):','Enter subject number (e.g. 01:'};
dlgtitle        = 'Info window';
dims            = [1 30];
definput        = {'beads','01'}; % this is a default input (this should change)
answer          = inputdlg(prompt,dlgtitle,dims,definput);

startup.answer  = answer;    % participant number and task name
getpath         = answer{1}; % usefull to read the correct task name and run the main script

switch getpath
    
    case 'beads'
        
        taskpath = fullfile(workingpath, 'beads');
        addpath(genpath(taskpath));
        main_beads % run main script
        
    case 'economic'
        
        taskpath = fullfile(workingpath, 'economic');
        addpath(genpath(taskpath));
        main_eco % run main script
        
    case 'face'
        
        taskpath = fullfile(workingpath, 'face');
        addpath(genpath(taskpath));
        main_face % run main script
    
end

clear startpath workingpath