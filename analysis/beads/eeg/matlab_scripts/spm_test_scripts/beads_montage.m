% beads montage creation

% load the data as SPM object
D                   = spm_eeg_load('spmeeg_sub_01_beads_block_01.mat');

montage.labelorg    = D.chanlabels; % define old label organisation

% the new label organisation should include EEG channels, a vertical and a horizontal EOG
montage.labelnew    = [montage.labelorg(1:64), 'HEOG', 'VEOG']; % new 
tra                 = eye(D.nchannels);                         % def identity matrix
tra(65:end, :)      = [];                                       % remove EXG rows 
tra                 = detrend(tra, 'constant');                 % use "detrend to average across sensors

tra(65, [67 68])    = [1 -1];                                   % HEOG
tra(66, [69 70])    = [1 -1];                                   % VEOG
montage.tra         = tra;

% save the montage 
save beads_montage.mat montage
