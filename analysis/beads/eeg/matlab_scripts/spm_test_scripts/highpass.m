spm('defaults', 'eeg');

%% define dataset S and run high pass filter 
S = [];
S.D = '/Users/christinadelta/Desktop/os_data/beads/spmdir/output/Mspmeeg_sub_01_beads_block_01.mat';
S.type = 'butterworth';
S.band = 'high';
S.freq = 0.1;
S.dir = 'twopass';
S.order = 5;
S.prefix = 'f';
D = spm_eeg_filter(S);


