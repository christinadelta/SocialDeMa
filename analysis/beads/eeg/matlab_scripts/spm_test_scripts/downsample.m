spm('defaults', 'eeg');

%% create dataset S and run downsample function 
S = [];
S.D = '/Users/christinadelta/Desktop/os_data/beads/spmdir/output/fMspmeeg_sub_01_beads_block_01.mat';
S.fsample_new = 256;
S.method = 'resample';
S.prefix = 'd';
D = spm_eeg_downsample(S);


