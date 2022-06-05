spm('defaults', 'eeg');

%% baseline correction
S = [];
S.D = '/Users/christinadelta/Desktop/os_data/beads/spmdir/output/efdfMspmeeg_sub_01_beads_block_01.mat';
S.timewin = [-500 -50];
S.prefix = 'b';
S.save = 1;
S.updatehistory = 1;
D = spm_eeg_bc(S);


