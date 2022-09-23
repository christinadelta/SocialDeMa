% % subjected specific MEEG objects
subout          = fullfile(outDir, sprintf('sub-%02d', sub));

S                       = [];
S.D                     = fullfile(subout, sprintf('maceerpfdfMspmeeg_sub_%02d_beads_block_01.mat', sub));
S.object                = spm_eeg_load(S.D);
tmp                     = S.object.conditions

% STATEMENT 1: if first cond is draw (either easy or difficult)
if size(tmp{1},2) == 8 
    
    % STATEMENT 2: if first condition is "draw"
    if tmp{1} == 'easydraw' % if first condition is easy
        
        if size(tmp{3},2) == 8 % if 3rd condition is diffdraw
            % [easy-draw easy-urn diff-draw diff-urn]
            % 1. Urn vs Draw 
            S                       = [];
            S.D                     = fullfile(subout, sprintf('maceerpfdfMspmeeg_sub_%02d_beads_block_01.mat', sub));
            S.c                     = [-1 1 -1 1];
            S.label                 = {'urnVSdraw'};
            S.weighted              = 1;
            S.prefix                = 'wud_';
            D                       = spm_eeg_contrast(S);

            % 2. Difficult vs Easy
            S                       = [];
            S.D                     = fullfile(subout, sprintf('maceerpfdfMspmeeg_sub_%02d_beads_block_01.mat', sub));
            S.c                     = [-1 -1 1 1];
            S.label                 = {'DiffVsEasy'};
            S.weighted              = 1;
            S.prefix                = 'wde_';
            D                       = spm_eeg_contrast(S);

            % 3. Interaction 
            S                       = [];
            S.D                     = fullfile(subout, sprintf('maceerpfdfMspmeeg_sub_%02d_beads_block_01.mat', sub));
            S.c                     = [-1 1 1 -1];
            S.label                 = {'interaction'};
            S.weighted              = 1;
            S.prefix                = 'wi_';
            D                       = spm_eeg_contrast(S);

            % 4. Only urns contrast
            S                       = [];
            S.D                     = fullfile(subout, sprintf('maceerpfdfMspmeeg_sub_%02d_beads_block_01.mat', sub));
            S.c                     = [0 1 0 1];
            S.label                 = {'onlyurn'};
            S.weighted              = 1;
            S.prefix                = 'wu_';
            D                       = spm_eeg_contrast(S);

            % 5. Only draws contrast
            S                       = [];
            S.D                     = fullfile(subout, sprintf('maceerpfdfMspmeeg_sub_%02d_beads_block_01.mat', sub));
            S.c                     = [1 0 1 0];
            S.label                 = {'onlydraw'};
            S.weighted              = 1;
            S.prefix                = 'wd_';
            D                       = spm_eeg_contrast(S);
            
        elseif size(tmp{3},2) == 7 % if 3rd condition is diffurn
            
            % [easy-draw easy-urn diff-urn diff-draw]
            % 1. Urn vs Draw 
            S                       = [];
            S.D                     = fullfile(subout, sprintf('maceerpfdfMspmeeg_sub_%02d_beads_block_01.mat', sub));
            S.c                     = [-1 1 1 -1];
            S.label                 = {'urnVSdraw'};
            S.weighted              = 1;
            S.prefix                = 'wud_';
            D                       = spm_eeg_contrast(S);

            % 2. Difficult vs Easy
            S                       = [];
            S.D                     = fullfile(subout, sprintf('maceerpfdfMspmeeg_sub_%02d_beads_block_01.mat', sub));
            S.c                     = [-1 -1 1 1];
            S.label                 = {'DiffVsEasy'};
            S.weighted              = 1;
            S.prefix                = 'wde_';
            D                       = spm_eeg_contrast(S);

            % 3. Interaction 
            S                       = [];
            S.D                     = fullfile(subout, sprintf('maceerpfdfMspmeeg_sub_%02d_beads_block_01.mat', sub));
            S.c                     = [-1 1 -1 1];
            S.label                 = {'interaction'};
            S.weighted              = 1;
            S.prefix                = 'wi_';
            D                       = spm_eeg_contrast(S);

            % 4. Only urns contrast
            S                       = [];
            S.D                     = fullfile(subout, sprintf('maceerpfdfMspmeeg_sub_%02d_beads_block_01.mat', sub));
            S.c                     = [0 1 1 0];
            S.label                 = {'onlyurn'};
            S.weighted              = 1;
            S.prefix                = 'wu_';
            D                       = spm_eeg_contrast(S);

            % 5. Only draws contrast
            S                       = [];
            S.D                     = fullfile(subout, sprintf('maceerpfdfMspmeeg_sub_%02d_beads_block_01.mat', sub));
            S.c                     = [1 0 0 1];
            S.label                 = {'onlydraw'};
            S.weighted              = 1;
            S.prefix                = 'wd_';
            D                       = spm_eeg_contrast(S);
 
        end
        
    % STAMENT 2: if first condition is draw
    elseif tmp{1} == 'diffdraw' % if first condition is difficult  
        
        if size(tmp{3},2) == 8 % if 3rd condition is easydraw
            
            % [diff-draw diff-urn easy-draw easy-urn]
            % 1. Urn vs Draw 
            S                       = [];
            S.D                     = fullfile(subout, sprintf('maceerpfdfMspmeeg_sub_%02d_beads_block_01.mat', sub));
            S.c                     = [-1 1 -1 1];
            S.label                 = {'urnVSdraw'};
            S.weighted              = 1;
            S.prefix                = 'wud_';
            D                       = spm_eeg_contrast(S);

            % 2. Difficult vs Easy
            S                       = [];
            S.D                     = fullfile(subout, sprintf('maceerpfdfMspmeeg_sub_%02d_beads_block_01.mat', sub));
            S.c                     = [1 1 -1 -1];
            S.label                 = {'DiffVsEasy'};
            S.weighted              = 1;
            S.prefix                = 'wde_';
            D                       = spm_eeg_contrast(S);

            % 3. Interaction 
            S                       = [];
            S.D                     = fullfile(subout, sprintf('maceerpfdfMspmeeg_sub_%02d_beads_block_01.mat', sub));
            S.c                     = [1 -1 -1 1];
            S.label                 = {'interaction'};
            S.weighted              = 1;
            S.prefix                = 'wi_';
            D                       = spm_eeg_contrast(S);

            % 4. Only urns contrast
            S                       = [];
            S.D                     = fullfile(subout, sprintf('maceerpfdfMspmeeg_sub_%02d_beads_block_01.mat', sub));
            S.c                     = [0 1 0 1];
            S.label                 = {'onlyurn'};
            S.weighted              = 1;
            S.prefix                = 'wu_';
            D                       = spm_eeg_contrast(S);

            % 5. Only draws contrast
            S                       = [];
            S.D                     = fullfile(subout, sprintf('maceerpfdfMspmeeg_sub_%02d_beads_block_01.mat', sub));
            S.c                     = [1 0 1 0];
            S.label                 = {'onlydraw'};
            S.weighted              = 1;
            S.prefix                = 'wd_';
            D                       = spm_eeg_contrast(S);
            
        elseif size(tmp{3},2) == 7 % if 3rd condition is easyurn
            
            % [diff-draw diff-urn easy-urn easy-draw]
            % 1. Urn vs Draw 
            S                       = [];
            S.D                     = fullfile(subout, sprintf('maceerpfdfMspmeeg_sub_%02d_beads_block_01.mat', sub));
            S.c                     = [-1 1 1 -1];
            S.label                 = {'urnVSdraw'};
            S.weighted              = 1;
            S.prefix                = 'wud_';
            D                       = spm_eeg_contrast(S);

            % 2. Difficult vs Easy
            S                       = [];
            S.D                     = fullfile(subout, sprintf('maceerpfdfMspmeeg_sub_%02d_beads_block_01.mat', sub));
            S.c                     = [1 1 -1 -1];
            S.label                 = {'DiffVsEasy'};
            S.weighted              = 1;
            S.prefix                = 'wde_';
            D                       = spm_eeg_contrast(S);

            % 3. Interaction 
            S                       = [];
            S.D                     = fullfile(subout, sprintf('maceerpfdfMspmeeg_sub_%02d_beads_block_01.mat', sub));
            S.c                     = [1 -1 1 -1];
            S.label                 = {'interaction'};
            S.weighted              = 1;
            S.prefix                = 'wi_';
            D                       = spm_eeg_contrast(S);

            % 4. Only urns contrast
            S                       = [];
            S.D                     = fullfile(subout, sprintf('maceerpfdfMspmeeg_sub_%02d_beads_block_01.mat', sub));
            S.c                     = [0 1 1 0];
            S.label                 = {'onlyurn'};
            S.weighted              = 1;
            S.prefix                = 'wu_';
            D                       = spm_eeg_contrast(S);

            % 5. Only draws contrast
            S                       = [];
            S.D                     = fullfile(subout, sprintf('maceerpfdfMspmeeg_sub_%02d_beads_block_01.mat', sub));
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
            
            % [easy-urn easy-draw diff-draw diff-urn]
            % 1. Urn vs Draw 
            S                       = [];
            S.D                     = fullfile(subout, sprintf('maceerpfdfMspmeeg_sub_%02d_beads_block_01.mat', sub));
            S.c                     = [1 -1 -1 1];
            S.label                 = {'urnVSdraw'};
            S.weighted              = 1;
            S.prefix                = 'wud_';
            D                       = spm_eeg_contrast(S);

            % 2. Difficult vs Easy
            S                       = [];
            S.D                     = fullfile(subout, sprintf('maceerpfdfMspmeeg_sub_%02d_beads_block_01.mat', sub));
            S.c                     = [-1 -1 1 1];
            S.label                 = {'DiffVsEasy'};
            S.weighted              = 1;
            S.prefix                = 'wde_';
            D                       = spm_eeg_contrast(S);

            % 3. Interaction 
            S                       = [];
            S.D                     = fullfile(subout, sprintf('maceerpfdfMspmeeg_sub_%02d_beads_block_01.mat', sub));
            S.c                     = [1 -1 1 -1];
            S.label                 = {'interaction'};
            S.weighted              = 1;
            S.prefix                = 'wi_';
            D                       = spm_eeg_contrast(S);

            % 4. Only urns contrast
            S                       = [];
            S.D                     = fullfile(subout, sprintf('maceerpfdfMspmeeg_sub_%02d_beads_block_01.mat', sub));
            S.c                     = [1 0 0 1];
            S.label                 = {'onlyurn'};
            S.weighted              = 1;
            S.prefix                = 'wu_';
            D                       = spm_eeg_contrast(S);

            % 5. Only draws contrast
            S                       = [];
            S.D                     = fullfile(subout, sprintf('maceerpfdfMspmeeg_sub_%02d_beads_block_01.mat', sub));
            S.c                     = [0 1 1 0];
            S.label                 = {'onlydraw'};
            S.weighted              = 1;
            S.prefix                = 'wd_';
            D                       = spm_eeg_contrast(S);
            
        elseif size(tmp{3},2) == 7 %
            
            % [easy-urn easy-draw diff-urn diff-draw]
            % 1. Urn vs Draw 
            S                       = [];
            S.D                     = fullfile(subout, sprintf('maceerpfdfMspmeeg_sub_%02d_beads_block_01.mat', sub));
            S.c                     = [1 -1 1 -1];
            S.label                 = {'urnVSdraw'};
            S.weighted              = 1;
            S.prefix                = 'wud_';
            D                       = spm_eeg_contrast(S);

            % 2. Difficult vs Easy
            S                       = [];
            S.D                     = fullfile(subout, sprintf('maceerpfdfMspmeeg_sub_%02d_beads_block_01.mat', sub));
            S.c                     = [-1 -1 1 1];
            S.label                 = {'DiffVsEasy'};
            S.weighted              = 1;
            S.prefix                = 'wde_';
            D                       = spm_eeg_contrast(S);

            % 3. Interaction 
            S                       = [];
            S.D                     = fullfile(subout, sprintf('maceerpfdfMspmeeg_sub_%02d_beads_block_01.mat', sub));
            S.c                     = [1 -1 -1 1];
            S.label                 = {'interaction'};
            S.weighted              = 1;
            S.prefix                = 'wi_';
            D                       = spm_eeg_contrast(S);

            % 4. Only urns contrast
            S                       = [];
            S.D                     = fullfile(subout, sprintf('maceerpfdfMspmeeg_sub_%02d_beads_block_01.mat', sub));
            S.c                     = [1 0 1 0];
            S.label                 = {'onlyurn'};
            S.weighted              = 1;
            S.prefix                = 'wu_';
            D                       = spm_eeg_contrast(S);

            % 5. Only draws contrast
            S                       = [];
            S.D                     = fullfile(subout, sprintf('maceerpfdfMspmeeg_sub_%02d_beads_block_01.mat', sub));
            S.c                     = [0 1 0 1];
            S.label                 = {'onlydraw'};
            S.weighted              = 1;
            S.prefix                = 'wd_';
            D                       = spm_eeg_contrast(S);

        end
        
    % STATEMENT 2: if first condition is "diff"    
    elseif tmp{1} == 'diffurn' % if first condition is diff urn
        
        if size(tmp{3},2) == 8 %
            
            % [diff-urn diff-draw easy-draw easy-urn]
            % 1. Urn vs Draw 
            S                       = [];
            S.D                     = fullfile(subout, sprintf('maceerpfdfMspmeeg_sub_%02d_beads_block_01.mat', sub));
            S.c                     = [1 -1 -1 1];
            S.label                 = {'urnVSdraw'};
            S.weighted              = 1;
            S.prefix                = 'wud_';
            D                       = spm_eeg_contrast(S);

            % 2. Difficult vs Easy
            S                       = [];
            S.D                     = fullfile(subout, sprintf('maceerpfdfMspmeeg_sub_%02d_beads_block_01.mat', sub));
            S.c                     = [1 1 -1 -1];
            S.label                 = {'DiffVsEasy'};
            S.weighted              = 1;
            S.prefix                = 'wde_';
            D                       = spm_eeg_contrast(S);

            % 3. Interaction 
            S                       = [];
            S.D                     = fullfile(subout, sprintf('maceerpfdfMspmeeg_sub_%02d_beads_block_01.mat', sub));
            S.c                     = [-1 1 -1 1];
            S.label                 = {'interaction'};
            S.weighted              = 1;
            S.prefix                = 'wi_';
            D                       = spm_eeg_contrast(S);

            % 4. Only urns contrast
            S                       = [];
            S.D                     = fullfile(subout, sprintf('maceerpfdfMspmeeg_sub_%02d_beads_block_01.mat', sub));
            S.c                     = [1 0 0 1];
            S.label                 = {'onlyurn'};
            S.weighted              = 1;
            S.prefix                = 'wu_';
            D                       = spm_eeg_contrast(S);

            % 5. Only draws contrast
            S                       = [];
            S.D                     = fullfile(subout, sprintf('maceerpfdfMspmeeg_sub_%02d_beads_block_01.mat', sub));
            S.c                     = [0 1 1 0];
            S.label                 = {'onlydraw'};
            S.weighted              = 1;
            S.prefix                = 'wd_';
            D                       = spm_eeg_contrast(S);
            
        elseif size(tmp{3},2) == 7 %
            
            % [diff-urn diff-draw easy-urn easy-draw]
            % 1. Urn vs Draw 
            S                       = [];
            S.D                     = fullfile(subout, sprintf('maceerpfdfMspmeeg_sub_%02d_beads_block_01.mat', sub));
            S.c                     = [1 -1 1 -1];
            S.label                 = {'urnVSdraw'};
            S.weighted              = 1;
            S.prefix                = 'wud_';
            D                       = spm_eeg_contrast(S);

            % 2. Difficult vs Easy
            S                       = [];
            S.D                     = fullfile(subout, sprintf('maceerpfdfMspmeeg_sub_%02d_beads_block_01.mat', sub));
            S.c                     = [1 1 -1 -1];
            S.label                 = {'DiffVsEasy'};
            S.weighted              = 1;
            S.prefix                = 'wde_';
            D                       = spm_eeg_contrast(S);

            % 3. Interaction 
            S                       = [];
            S.D                     = fullfile(subout, sprintf('maceerpfdfMspmeeg_sub_%02d_beads_block_01.mat', sub));
            S.c                     = [-1 1 1 -1];
            S.label                 = {'interaction'};
            S.weighted              = 1;
            S.prefix                = 'wi_';
            D                       = spm_eeg_contrast(S);

            % 4. Only urns contrast
            S                       = [];
            S.D                     = fullfile(subout, sprintf('maceerpfdfMspmeeg_sub_%02d_beads_block_01.mat', sub));
            S.c                     = [1 0 1 0];
            S.label                 = {'onlyurn'};
            S.weighted              = 1;
            S.prefix                = 'wu_';
            D                       = spm_eeg_contrast(S);

            % 5. Only draws contrast
            S                       = [];
            S.D                     = fullfile(subout, sprintf('maceerpfdfMspmeeg_sub_%02d_beads_block_01.mat', sub));
            S.c                     = [0 1 0 1];
            S.label                 = {'onlydraw'};
            S.weighted              = 1;
            S.prefix                = 'wd_';
            D                       = spm_eeg_contrast(S);
   
        end

    end % end of main statement 2
 
end % end of main if statement