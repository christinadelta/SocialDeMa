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
% Beads task. These are:

% 1. createMontage.m -- This function creates a montage (performs channel
% selection) and re-refernces the signal to the average of all electrodes.
% In pre-processing we reference the data to the average of all electrodes.
% This is one of the most common and well known methods. This method
% requires excluding noisy electrodes from referencing. Use this log file
% to see if a given participant has a bad channel. In case there is a bad
% channel exclude it from the createMontage.m function (for more info see
% the documentation within the createMontage.m file).

% 2. beads_trialdef.m -- this function re-writes the events as as draw
% choices and urn choices (for the 0.8 & 0.6 conditions). During the
% recording there was no way to code the triggers as draw and urn choices
% (because the trigger is sent right after the bead presentation screen is
% flipped), thus we need to do this change "manually" before epoching. For
% more info on how this is accomplished, see the documentation of the 
% beads_trialdef.m file. 

 


