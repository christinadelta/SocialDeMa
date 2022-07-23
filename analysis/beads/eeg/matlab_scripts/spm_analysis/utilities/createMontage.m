function S = createMontage(S)

% this function is part of the beads EEG preprocessing/analysis script. 
% Input: S struct with the MEEG object from the conversion step 
% Output: updated S struct including path for the montage .mat file that
% this func creates 

% this function requires knowing in advance if the current subject dataset
% has a noisy electrode (e.g. pilot sub has one noisy channel [channel
% 25]). This will be excluded.

% TODO:
% 1. re-referencing needs fixing

%%% -------------------- Run the function
% load the object 
D                   = S.D;
obj                 = spm_eeg_load(D);

montage.labelorg    = obj.chanlabels; % old labels
montage.labelnew    = [montage.labelorg(1:64), 'HEOG', 'VEOG']; % new labels

% re-reference 
tra                 = eye(obj.nchannels);
tra(65:end, :)      = []; % remove the last two EXG channels (EXG7, EXG8)

% also exclude channel 25 (PO7) only for this subject. If a subject doesn't have a noisy channel, comment this line
% tra(25, :)          = []; 
tra                 = detrend(tra, 'constant'); % re-reference

% HEOG
tra(64, [67 68])    = [1 -1];

% VEOG
tra(65, [69 70])    = [1 -1];

montage.tra         = tra;

% save montage in jobs dir
S.sjob              = fullfile(S.jobpath, 'montage.mat');
save(S.sjob, 'montage');

return