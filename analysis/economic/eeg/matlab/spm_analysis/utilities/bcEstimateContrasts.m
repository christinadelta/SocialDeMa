function [] = bcEstimateContrasts(S,tmp,sub)

% created on 15/01/2023 for economic (and facial-attractiveness) task preprocessing
% the function runs from the main preprocessing/analysis script
% Inputs:
%           - S structure (with filename)
%           - tmp (cell that contains the order of the conditions for that subject)
%           - subject number (this is needed for saving the new file in the correct output directory)
%           

% Saves the new meeg object (contrast) in the subject's output dir
% extract filename from S struct
fn = S.D;

% -----------------------------
% extract filename from S struct
S                       = [];
S.D                     = fn;
S.c                     = [-1 1];
S.label                 = {'rejectVSaccept'};
S.weighted              = 1;
S.prefix                = 'wra_';
D                       = spm_eeg_contrast(S);


end 