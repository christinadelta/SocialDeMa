function [] = beadsEstimateContrast(S,tmp,sub)

% created on 23/09/2022 for beads task preprocessing
% the function runs from the main preprocessing/analysis script
% Inputs:
%           - S structure (with filename)
%           - tmp (cell that contains the order of the conditions for that subject)
%           - subject number (this is needed for saving the new file in the correct output directory)
%           

% Saves the new meeg object (contrast) in the subject's output dir

% extract filename from S struct
fn = S.D;

%--------------------------------------------------------------

% STATEMENT 1: if first cond is draw (either easy or difficult)
if size(tmp{1},2) == 8 
    
    % STATEMENT 2: if first condition is "easy"
    if tmp{1} == 'easydraw' % if first condition is easydraw
        
        if size(tmp{3},2) == 8 % if 3rd condition is diffdraw
            
            % order 1: [easy-draw easy-urn diff-draw diff-urn]
            fprintf('order 1') % just for diagnostics
            
            % 1. Urn vs Draw 
            S                       = [];
            S.D                     = fn;
            S.c                     = [-1 1 -1 1];
            S.label                 = {'urnVSdraw'};
            S.weighted              = 1;
            S.prefix                = 'wud_';
            D                       = spm_eeg_contrast(S);

            % 2. Difficult vs Easy
            S                       = [];
            S.D                     = fn;
            S.c                     = [-1 -1 1 1];
            S.label                 = {'DiffVsEasy'};
            S.weighted              = 1;
            S.prefix                = 'wde_';
            D                       = spm_eeg_contrast(S);

            % 3. Interaction 
            S                       = [];
            S.D                     = fn;
            S.c                     = [-1 1 1 -1];
            S.label                 = {'interaction'};
            S.weighted              = 1;
            S.prefix                = 'wi_';
            D                       = spm_eeg_contrast(S);

            % 4. Only urns contrast
            S                       = [];
            S.D                     = fn;
            S.c                     = [0 1 0 1];
            S.label                 = {'onlyurn'};
            S.weighted              = 1;
            S.prefix                = 'wu_';
            D                       = spm_eeg_contrast(S);

            % 5. Only draws contrast
            S                       = [];
            S.D                     = fn;
            S.c                     = [1 0 1 0];
            S.label                 = {'onlydraw'};
            S.weighted              = 1;
            S.prefix                = 'wd_';
            D                       = spm_eeg_contrast(S);
            
        elseif size(tmp{3},2) == 7 % if 3rd condition is diffurn
            
            % order 2: [easy-draw easy-urn diff-urn diff-draw]
            fprintf('order 2')
            
            % 1. Urn vs Draw 
            S                       = [];
            S.D                     = fn;
            S.c                     = [-1 1 1 -1];
            S.label                 = {'urnVSdraw'};
            S.weighted              = 1;
            S.prefix                = 'wud_';
            D                       = spm_eeg_contrast(S);

            % 2. Difficult vs Easy
            S                       = [];
            S.D                     = fn;
            S.c                     = [-1 -1 1 1];
            S.label                 = {'DiffVsEasy'};
            S.weighted              = 1;
            S.prefix                = 'wde_';
            D                       = spm_eeg_contrast(S);

            % 3. Interaction 
            S                       = [];
            S.D                     = fn;
            S.c                     = [-1 1 -1 1];
            S.label                 = {'interaction'};
            S.weighted              = 1;
            S.prefix                = 'wi_';
            D                       = spm_eeg_contrast(S);

            % 4. Only urns contrast
            S                       = [];
            S.D                     = fn;
            S.c                     = [0 1 1 0];
            S.label                 = {'onlyurn'};
            S.weighted              = 1;
            S.prefix                = 'wu_';
            D                       = spm_eeg_contrast(S);

            % 5. Only draws contrast
            S                       = [];
            S.D                     = fn;
            S.c                     = [1 0 0 1];
            S.label                 = {'onlydraw'};
            S.weighted              = 1;
            S.prefix                = 'wd_';
            D                       = spm_eeg_contrast(S);
 
        end
        
    % STAMENT 2: if first condition is diff
    elseif tmp{1} == 'diffdraw' % if first condition is diffdraw 
        
        if size(tmp{3},2) == 8 % if 3rd condition is easydraw
            
            % order 3; [diff-draw diff-urn easy-draw easy-urn]
            fprintf('order 3')
            
            % 1. Urn vs Draw 
            S                       = [];
            S.D                     = fn;
            S.c                     = [-1 1 -1 1];
            S.label                 = {'urnVSdraw'};
            S.weighted              = 1;
            S.prefix                = 'wud_';
            D                       = spm_eeg_contrast(S);

            % 2. Difficult vs Easy
            S                       = [];
            S.D                     = fn;
            S.c                     = [1 1 -1 -1];
            S.label                 = {'DiffVsEasy'};
            S.weighted              = 1;
            S.prefix                = 'wde_';
            D                       = spm_eeg_contrast(S);

            % 3. Interaction 
            S                       = [];
            S.D                     = fn;
            S.c                     = [1 -1 -1 1];
            S.label                 = {'interaction'};
            S.weighted              = 1;
            S.prefix                = 'wi_';
            D                       = spm_eeg_contrast(S);

            % 4. Only urns contrast
            S                       = [];
            S.D                     = fn;
            S.c                     = [0 1 0 1];
            S.label                 = {'onlyurn'};
            S.weighted              = 1;
            S.prefix                = 'wu_';
            D                       = spm_eeg_contrast(S);

            % 5. Only draws contrast
            S                       = [];
            S.D                     = fn;
            S.c                     = [1 0 1 0];
            S.label                 = {'onlydraw'};
            S.weighted              = 1;
            S.prefix                = 'wd_';
            D                       = spm_eeg_contrast(S);
            
        elseif size(tmp{3},2) == 7 % if 3rd condition is easyurn
            
            % order 4: [diff-draw diff-urn easy-urn easy-draw]
            fprintf('order 4')
            
            % 1. Urn vs Draw 
            S                       = [];
            S.D                     = fn;
            S.c                     = [-1 1 1 -1];
            S.label                 = {'urnVSdraw'};
            S.weighted              = 1;
            S.prefix                = 'wud_';
            D                       = spm_eeg_contrast(S);

            % 2. Difficult vs Easy
            S                       = [];
            S.D                     = fn;
            S.c                     = [1 1 -1 -1];
            S.label                 = {'DiffVsEasy'};
            S.weighted              = 1;
            S.prefix                = 'wde_';
            D                       = spm_eeg_contrast(S);

            % 3. Interaction 
            S                       = [];
            S.D                     = fn;
            S.c                     = [1 -1 1 -1];
            S.label                 = {'interaction'};
            S.weighted              = 1;
            S.prefix                = 'wi_';
            D                       = spm_eeg_contrast(S);

            % 4. Only urns contrast
            S                       = [];
            S.D                     = fn;
            S.c                     = [0 1 1 0];
            S.label                 = {'onlyurn'};
            S.weighted              = 1;
            S.prefix                = 'wu_';
            D                       = spm_eeg_contrast(S);

            % 5. Only draws contrast
            S                       = [];
            S.D                     = fn;
            S.c                     = [1 0 0 1];
            S.label                 = {'onlydraw'};
            S.weighted              = 1;
            S.prefix                = 'wd_';
            D                       = spm_eeg_contrast(S);
            
        end
    end % end of main statemrnt 2

% STATEMENT 1: if first cond is urn (either easy or difficult)
elseif size(tmp{1},2) == 7
    
    % STATEMENT 2: if first condition is "easy"
    if tmp{1} == 'easyurn' % if first condition is easy urn
        if size(tmp{3},2) == 8 %
            
            % order 5: [easy-urn easy-draw diff-draw diff-urn]
            fprintf('order 5')
            
            % 1. Urn vs Draw 
            S                       = [];
            S.D                     = fn;
            S.c                     = [1 -1 -1 1];
            S.label                 = {'urnVSdraw'};
            S.weighted              = 1;
            S.prefix                = 'wud_';
            D                       = spm_eeg_contrast(S);

            % 2. Difficult vs Easy
            S                       = [];
            S.D                     = fn;
            S.c                     = [-1 -1 1 1];
            S.label                 = {'DiffVsEasy'};
            S.weighted              = 1;
            S.prefix                = 'wde_';
            D                       = spm_eeg_contrast(S);

            % 3. Interaction 
            S                       = [];
            S.D                     = fn;
            S.c                     = [1 -1 1 -1];
            S.label                 = {'interaction'};
            S.weighted              = 1;
            S.prefix                = 'wi_';
            D                       = spm_eeg_contrast(S);

            % 4. Only urns contrast
            S                       = [];
            S.D                     = fn;
            S.c                     = [1 0 0 1];
            S.label                 = {'onlyurn'};
            S.weighted              = 1;
            S.prefix                = 'wu_';
            D                       = spm_eeg_contrast(S);

            % 5. Only draws contrast
            S                       = [];
            S.D                     = fn;
            S.c                     = [0 1 1 0];
            S.label                 = {'onlydraw'};
            S.weighted              = 1;
            S.prefix                = 'wd_';
            D                       = spm_eeg_contrast(S);
            
        elseif size(tmp{3},2) == 7 %
            
            % order 6: [easy-urn easy-draw diff-urn diff-draw]
            fprintf('order 6')
            
            % 1. Urn vs Draw 
            S                       = [];
            S.D                     = fn;
            S.c                     = [1 -1 1 -1];
            S.label                 = {'urnVSdraw'};
            S.weighted              = 1;
            S.prefix                = 'wud_';
            D                       = spm_eeg_contrast(S);

            % 2. Difficult vs Easy
            S                       = [];
            S.D                     = fn;
            S.c                     = [-1 -1 1 1];
            S.label                 = {'DiffVsEasy'};
            S.weighted              = 1;
            S.prefix                = 'wde_';
            D                       = spm_eeg_contrast(S);

            % 3. Interaction 
            S                       = [];
            S.D                     = fn;
            S.c                     = [1 -1 -1 1];
            S.label                 = {'interaction'};
            S.weighted              = 1;
            S.prefix                = 'wi_';
            D                       = spm_eeg_contrast(S);

            % 4. Only urns contrast
            S                       = [];
            S.D                     = fn;
            S.c                     = [1 0 1 0];
            S.label                 = {'onlyurn'};
            S.weighted              = 1;
            S.prefix                = 'wu_';
            D                       = spm_eeg_contrast(S);

            % 5. Only draws contrast
            S                       = [];
            S.D                     = fn;
            S.c                     = [0 1 0 1];
            S.label                 = {'onlydraw'};
            S.weighted              = 1;
            S.prefix                = 'wd_';
            D                       = spm_eeg_contrast(S);

        end
        
    % STATEMENT 2: if first condition is "diff"    
    elseif tmp{1} == 'diffurn' % if first condition is diff urn
        
        if size(tmp{3},2) == 8 %
            
            % order 7: [diff-urn diff-draw easy-draw easy-urn]
            fprintf('order 7')
            
            % 1. Urn vs Draw 
            S                       = [];
            S.D                     = fn;
            S.c                     = [1 -1 -1 1];
            S.label                 = {'urnVSdraw'};
            S.weighted              = 1;
            S.prefix                = 'wud_';
            D                       = spm_eeg_contrast(S);

            % 2. Difficult vs Easy
            S                       = [];
            S.D                     = fn;
            S.c                     = [1 1 -1 -1];
            S.label                 = {'DiffVsEasy'};
            S.weighted              = 1;
            S.prefix                = 'wde_';
            D                       = spm_eeg_contrast(S);

            % 3. Interaction 
            S                       = [];
            S.D                     = fn;
            S.c                     = [-1 1 -1 1];
            S.label                 = {'interaction'};
            S.weighted              = 1;
            S.prefix                = 'wi_';
            D                       = spm_eeg_contrast(S);

            % 4. Only urns contrast
            S                       = [];
            S.D                     = fn;
            S.c                     = [1 0 0 1];
            S.label                 = {'onlyurn'};
            S.weighted              = 1;
            S.prefix                = 'wu_';
            D                       = spm_eeg_contrast(S);

            % 5. Only draws contrast
            S                       = [];
            S.D                     = fn;
            S.c                     = [0 1 1 0];
            S.label                 = {'onlydraw'};
            S.weighted              = 1;
            S.prefix                = 'wd_';
            D                       = spm_eeg_contrast(S);
            
        elseif size(tmp{3},2) == 7 %
            
            % order 8: [diff-urn diff-draw easy-urn easy-draw]
            fprintf('order 8')
            
            % 1. Urn vs Draw 
            S                       = [];
            S.D                     = fn;
            S.c                     = [1 -1 1 -1];
            S.label                 = {'urnVSdraw'};
            S.weighted              = 1;
            S.prefix                = 'wud_';
            D                       = spm_eeg_contrast(S);

            % 2. Difficult vs Easy
            S                       = [];
            S.D                     = fn;
            S.c                     = [1 1 -1 -1];
            S.label                 = {'DiffVsEasy'};
            S.weighted              = 1;
            S.prefix                = 'wde_';
            D                       = spm_eeg_contrast(S);

            % 3. Interaction 
            S                       = [];
            S.D                     = fn;
            S.c                     = [-1 1 1 -1];
            S.label                 = {'interaction'};
            S.weighted              = 1;
            S.prefix                = 'wi_';
            D                       = spm_eeg_contrast(S);

            % 4. Only urns contrast
            S                       = [];
            S.D                     = fn;
            S.c                     = [1 0 1 0];
            S.label                 = {'onlyurn'};
            S.weighted              = 1;
            S.prefix                = 'wu_';
            D                       = spm_eeg_contrast(S);

            % 5. Only draws contrast
            S                       = [];
            S.D                     = fn;
            S.c                     = [0 1 0 1];
            S.label                 = {'onlydraw'};
            S.weighted              = 1;
            S.prefix                = 'wd_';
            D                       = spm_eeg_contrast(S);
   
        end

    end % end of main statement 2
 
end % end of main if statement

return