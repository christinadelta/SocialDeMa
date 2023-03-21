% % PRE-PROCESSING AND ANALYSIS SCRIPT FOR THE BEADS TASK VERSION 1

% Part of the Optimal Stopping Problems Project

% CREATED: 21/03/2023
% This version will take over the previous (beads_prepro_v6.m) one which now is
% archived. 

% changes introduced in this version: 
% 1) I run everything using this script. The steps are described below 
% 2) I fixed the beta thingy after discussing with Nick 

%% PREPROCESSING STEPS %%

% 1. extract sub block data [logs, sequences]
% 2. store sub data based on conditions
% 3. average sub draws and acc for each condition

%% ANALYSES STEPS %%

% 1. run ideal observer (IO)
% 2. extract IO's sumpling rate and performance and average 
% 3. run statistical analyses using uaing the runBehavStats.m file:
%       a) run 2x2 anovas with factors: [probabilities, agent]
%       b) run pairwise comparisons to look at differences in model-human, 0.8 & 0.6 conditions

% 4. fit models:
%   model 1: one free params [beta]
%   model 2: two free params [beta & Cs]
%   parameter recovery
%   model comparison

% 5. compute AQ differences [draw-again higher-option]
% 6. next analysis steps will be introduced soon...

%% PLOTING %%

% 1. plot human behaviour and IO
% 2. plot human, IO and model fit behviours 
% 3. plot model fitting/parameter recovery stuff
% 4. plot regression scatterplots? -- for later 


%% IMPORTANT NOTE %%

% I save two different types of mat files. In 4 mat files I save the block
% information. This kind of file contains the trials of the given block, the
% number of draws for each trial, the balance of this trial (reward/loss),
% response, accuracy, etc..

% The rest mat files contain sequence/trial information, such as: each draw
% of the given trial, trial start, bead onset, response time, etc...

% TOTAL MAT FILES for each subject: 56
% BLOCK MAT FILES: 4 logs - USED for preprocessing
% SEQUENCE MAT FILES: 52 logs - NOT USED

% END OF PREAMBLE 

%% 

