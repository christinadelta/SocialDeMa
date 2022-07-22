%% PREPROCESSING & ANALYSIS OF BEADS EEG DATA

% preprocessing/analysis script was created in July 2022.
% VERSION 1 of formal analysis plan of the beads EEG data with SPM12 

% Details of Preprocessing steps can be found in the PDF file:
% https://docs.google.com/document/d/1xLbgGW23Dk4S0rfzSJbUMTxc2_srCTFf77gUFijCR5Y/edit#heading=h.3ewqtkwgw38n
% 
% The EEG data for beads task were recorded in four .bdf files -- one file
% for each block.
% for preprocessing, up to epoching, each file is pre-processed seperately 

%%%%% VERY IMPORTANT %%%%%: 
% there were cases that a participant would have a noisy electrode. Quality
% of the EEG data during recording was logged in this file:
% https://docs.google.com/spreadsheets/d/16vhhOgA19vZDzn-K2K8PTeD3FnRPTRPAypzS450-NSk/edit#gid=0

% all preprocessing steps and analyses are done using this script. Meaning that we call all
% the SPM functions using this script.
% There shouldn't be a need to run any of the functions using the SMP eeg
% GUI. 

% The are two functions that are modified or created specifically for the
% Beads task in the 'utilities' dir. These are:

% 1. createMontage.m -- This function creates a montage (performs channel
% selection) and re-refernces the signal to the average of all electrodes.
% In pre-processing we reference the data to the average of all electrodes.
% This is one of the most common and well known methods. This method
% requires excluding noisy electrodes from referencing. Use this log file
% to see if a given participant has a bad channel. In case there is a bad
% channel exclude it from the createMontage.m function (for more info see
% the documentation within the createMontage.m file).

% 2. beadsTrialdef.m -- this function re-writes the events as as draw
% choices and urn choices (for the 0.8 & 0.6 conditions). During the
% recording there was no way to code the triggers as draw and urn choices
% (because the trigger is sent right after the bead presentation screen is
% flipped), thus we need to do this change "manually" before epoching. For
% more info on how this is accomplished, see the documentation of the 
% beadsTrialdef.m file. 

%% Create required directories and define paths %%

% Add SPM12 to the matlab path if needed:
% addpath /Users/christinadelta/neuro_software/spm12 % change this to your
% own SPM12 path
% savepath

% create spm directories that will be used for storing the output MEEG
% objects mat files and images. MEEG objects are stored in the output dir,
% .mat files in the jobs dir and images are stored in directories that are
% created automatically during conversion of the MEEG objects to .nii
% files.

basedir     = pwd;
spmdir      = fullfile(basedir, 'spmDir');
addpath(genpath(spmdir));

% if output and jobs directories do not exist, create them:
erpDir      = fullfile(spmdir, 'outerps');
tfrDir      = fullfile(spmdir, 'outtfrs');
jobsDir     = fullfile(spmdir, 'jobs');

if ~exist(erpDir, 'dir') && ~exist(tfrDir, 'dir') && ~exist(jobsDir, 'dir')
    mkdir(erpDir)
    mkdir(tfrDir)
    mkdir(jobsDir)
end

% set data path
datadir     = '/Users/christinadelta/Desktop/os_data/beads/subs/';  % change to your own data path 
subs        = dir(fullfile(datadir, '*sub*'));                      % how many subject folders?
nsubs       = length(subs);                                         % how many subjects?

taskname    = 'beads';
blocks      = 4;
conditions  = 2;
choices     = 2;

%% Preprocessing steps 1 - 8 %%

% these steps are used to preprocess the raw .bdf files, for each
% participant. Given that we will run evoked analysis and time-frequency
% analysis, we will low-pass filter the data in two ways, one for evoked
% analysis and one for TF analysis. Both MEEG objects will be preprocessed
% up until merging, however, only the ERPs MEEG object will undergo
% artefact rejection (the last step of preprocessing).

% loop over subjects 
for sub = 1:nsubs
    
    fprintf('loading beads block data\n')  
    subject = subs(sub).name;
    subdir  = fullfile(datadir,subject);
    fprintf('\t reading data from subject %d\n',sub); 
    
    % loop over blocks 
    for block = 1:blocks
        
        fprintf('\t\t loading block %d\n\n',block);
        
    end % end of blocks loop    
 
end % end of subjects loop

