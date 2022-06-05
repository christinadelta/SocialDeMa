spm('defaults', 'eeg');

%% create S dataset and run filter function set to lowpass 
S = [];
S.D = '/Users/christinadelta/Desktop/os_data/beads/spmdir/output/dfMspmeeg_sub_01_beads_block_01.mat';
S.type = 'butterworth';
S.band = 'low';
S.freq = 30;
S.dir = 'twopass';
S.order = 5;
S.prefix = 'f';
D = spm_eeg_filter(S);


