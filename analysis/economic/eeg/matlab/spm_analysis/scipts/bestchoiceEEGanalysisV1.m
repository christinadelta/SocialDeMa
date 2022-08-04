%% PREPROCESSING & ANALYSIS OF ECONOMIC AND FACIAL ATTRACTIVENESS (BEST-CHOICE) EEG DATA

% Preprocessing/analysis script was created in August 2022.
% VERSION 1 of formal analysis plan of the best-choice tasks EEG data with SPM12 

% BEST-CHOICE TASKS:
% 1. economic (best-choice) task
% 2. facial attractiveness (best-choice) task

% Details of Preprocessing steps can be found in the doc file:
% https://docs.google.com/document/d/1JAkRF-O4qCPNAXVtlYUEgGkUumYOu4-yJQ0xP9ycKbM/edit?usp=sharing
% 
% The EEG data for the best-choice tasks were recorded in two .bdf files -- one file
% for each block.
% for preprocessing, up to epoching, each file is pre-processed seperately 

%%%%% VERY IMPORTANT %%%%%: 
% there were cases that a participant would have a noisy electrode. Quality
% of the EEG data during recording was logged in this file:
% https://docs.google.com/spreadsheets/d/16vhhOgA19vZDzn-K2K8PTeD3FnRPTRPAypzS450-NSk/edit#gid=0

% all preprocessing steps and analyses are done using this script. Meaning that we call all
% the SPM functions using this script.
% There shouldn't be a need to run any of the functions using the SPM eeg
% GUI. 

% The are two functions that are modified or created specifically for the
% best-choice tasks in the 'utilities' dir. These are:

% 1. createMontage.m -- This function creates a montage (performs channel
% selection) and re-refernces the signal to the average of all electrodes.
% In pre-processing we reference the data to the average of all electrodes.
% This is one of the most common and well known methods. This method
% requires excluding noisy electrodes from referencing. Use this log file
% to see if a given participant has a bad channel. In case there is a bad
% channel exclude it from the createMontage.m function (for more info see
% the documentation within the createMontage.m file).

% 2. beadsTrialdef.m -- this function re-writes the events as as sample
% choices (reject choices) and accept choices. During the
% recording there was no way to code the triggers as reject and accept choices
% (because the trigger is sent right after the bead presentation screen is
% flipped), thus we need to do this change "manually" before epoching. For
% more info on how this is accomplished, see the documentation of the 
% bestchoiceTrialdef.m file. 

%%%% TODO:
% 1. FIX MONTAGE SPECIFICATION

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

% clear workspace
clear all
clc

basedir     = pwd;
% spmdir      = fullfile(basedir, 'spmDir');
% addpath(genpath(basedir));
addpath(genpath(fullfile(basedir, 'scripts')));
addpath(genpath(fullfile(basedir, 'utilities')));

% if output and jobs directories do not exist, create them:
outDir      = fullfile(basedir, 'output');
jobsDir     = fullfile(basedir, 'jobs');

% which task is it?
tasknum     = 1; % if 1=economic, if 2=facial attractivenes

if ~exist(outDir, 'dir') && ~exist(jobsDir, 'dir')
    mkdir(outDir)
    mkdir(jobsDir)
end

addpath(genpath(fullfile(basedir, 'jobs')));

if tasknum == 1
    taskname    = 'economic';
else
    taskname    = 'face';
end

datadir         = '/Users/christinadelta/Desktop/os_data/';             % change to your own path
taskdir         = fullfile(datadir, taskname, 'subs');

% set data path 
subs            = dir(fullfile(taskdir, '*sub*'));                      % how many subject folders?
nsubs           = length(subs);                                         % how many subjects?

blocks          = 2;
conditions      = 2;
choices         = 2;
analysest       = 2; % analyses types (erp & tf)
nconstrasts     = 3;
% contrastpref    = {'wud_' 'wde_' 'wi' 'wu_' 'wd_'};

%% Preprocessing steps [ERPs: 1 - 12, TFRs: 1 - 12] %%

% these steps are used to preprocess the raw .bdf files, for each
% participant. Given that we will run evoked analysis and time-frequency
% analysis, we will low-pass filter the data in two ways, one for evoked
% analysis and one for TF analysis. Both MEEG objects will be preprocessed
% up until merging, however, only the ERPs MEEG object will undergo
% artefact rejection (the last step of preprocessing).

% init spm 
spm('defaults', 'eeg');

% loop over subjects 
for sub = 1:nsubs
    
    fprintf('loading best-choice block data\n')  
    subject         = subs(sub).name;
    subdir          = fullfile(taskdir,subject);
    fprintf('\t reading data from subject %d\n',sub); 
    
    % create a subject sub-directory in outerps & outtfrs to store
    % subjected specific MEEG objects
    subout          = fullfile(outDir, sprintf('sub-%02d', sub));
    subjobs         = fullfile(jobsDir, sprintf('sub-%02d', sub));
    
    if ~exist(subout, 'dir') && ~exist(subjobs, 'dir')
        mkdir(subout)
        mkdir(subjobs)
    end
    
    % loop over blocks 
    for block = 1:blocks
        
        fprintf('\t\t loading block %d\n\n',block);
        blockfile       = fullfile(subdir, sprintf('sub_%02d_%s_block_%02d.bdf', sub, taskname, block));
        
        %% STEP 1. convert the bdf file to MEEG object
        % create S struct for conversion
        S               = [];
        S.dataset       = blockfile;
        S.mode          = 'continuous';
        S.channels      = {'Fp1', 'AF7', 'AF3', 'F1', 'F3', 'F5', 'F7', 'FT7', 'FC5', 'FC3', 'FC1', 'C1', 'C3', 'C5', 'T7',...
            'TP7', 'CP5', 'CP3', 'CP1', 'P1', 'P3', 'P5', 'P7', 'P9', 'PO7', 'PO3', 'O1', 'Iz', 'Oz', 'POz', 'Pz', 'CPz',...
            'Fpz', 'Fp2', 'AF8', 'AF4', 'AFz', 'Fz', 'F2', 'F4', 'F6', 'F8', 'FT8', 'FC6', 'FC4', 'FC2', 'FCz','Cz', 'C2',...
            'C4', 'C6', 'T8', 'TP8', 'CP6', 'CP4', 'CP2', 'P2', 'P4', 'P6', 'P8', 'P10', 'PO8',  'PO4', 'O2', 'EXG1',...
            'EXG2', 'EXG3', 'EXG4', 'EXG5', 'EXG6', 'EXG7', 'EXG8'};

        % convert bdf file to spm object and D struct
        S.eventpadding      = 0;
        % S.blocksize         = 3276800;
        S.checkboundary     = 1;
        S.saveorigheader    = 0;
        S.outpath           = fullfile(subout, sprintf('spmeeg_sub_%02d_%s_block_%02d.mat', sub, taskname, block));
        S.outfile           = S.outpath;
        S.timewin           = [];
        S.conditionlabels   = {'Undefined'};
        S.inputformat       = [];
        D                   = spm_eeg_convert(S); % convert raw data
        
        %% STEP 2. Create montage & re-reference 
        % create S struct for 
        S                   = [];
        S.D                 = fullfile(subout, sprintf('spmeeg_sub_%02d_%s_block_%02d.mat', sub, taskname, block));
        S.jobpath           = subjobs;
        S.block             = block;
        
        % Run createMontage.m function.
        % This function re-references by averaging across all electrodes. However,
        % an initial step requires knowing in advance if there is a noisy channel
        % and exclude it from averaging. Montage creation thus, need to be done
        % individually for every subject (e.g. pilot sub has one noisy channel [channel
        % 25]. This needs to be removed from averaging when re-referencing.
        % The createMontage.m file is in basedir/utilities; the function
        % needs to be modified manually 
        S                   = createMontageBestChoice(S);
        
        S.mode              = 'write';
        S.blocksize         = 655360;
        S.prefix            = 'M';
        S.montage           = fullfile(S.jobpath, 'montage.mat');
        S.keepothers        = 0;
        S.keepsensors       = 1;
        S.updatehistory     = 1;
        D                   = spm_eeg_montage(S);
        
        %% STEP 3. High-pass filter 
        % Init S struct
        S                   = [];
        S.D                 = fullfile(subout, sprintf('Mspmeeg_sub_%02d_%s_block_%02d.mat', sub, taskname, block));
        S.type              = 'butterworth';
        S.band              = 'high';
        S.freq              = 0.1;
        S.dir               = 'twopass';
        S.order             = 5;
        S.prefix            = 'f';
        D                   = spm_eeg_filter(S);
        
        %% STEP 4. Downsample
        % Init S struct
        S                   = [];
        S.D                 = fullfile(subout, sprintf('fMspmeeg_sub_%02d_%s_block_%02d.mat', sub, taskname, block));
        S.fsample_new       = 256;
        S.method            = 'resample';
        S.prefix            = 'd';
        D                   = spm_eeg_downsample(S);
        
        %% STEP 5. Low-pass filter 
        
        % at this stage we will low-pass filter the downsampled file twice:
        % 1. for ERP analysis with frequency cutoff: 30Hz
        % 2. for TF analysis with frequency cutoff: 110Hz
        
        for i = 1:analysest
            S               = [];
            S.D             = fullfile(subout, sprintf('dfMspmeeg_sub_%02d_%s_block_%02d.mat', sub, taskname, block));
            S.type          = 'butterworth';
            S.band          = 'low';
            S.dir           = 'twopass';
            S.order         = 5;
            
            if i == 1 % if it's for ERP analysis
                S.freq      = 30;
                S.prefix    = 'erpf';
                
            else % if it's tf analysis
                S.freq      = 110;
                S.prefix    = 'tfrf';
            end
            
            D               = spm_eeg_filter(S);
            
        end % end of analysis types loop
        
        
        
    
    end % end of blocks loop
    
end % end of subjects loop




