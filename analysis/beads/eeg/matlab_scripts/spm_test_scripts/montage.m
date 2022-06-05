spm('defaults', 'eeg');

%% run the "beads_montage.m" script/function that creates montage.mat

%% define dataset S and run montage from preprocessing 
S = [];
S.D = '/Users/christinadelta/Desktop/os_data/beads/spmdir/output/spmeeg_sub_01_beads_block_01.mat';
S.mode = 'write';
S.blocksize = 655360;
S.prefix = 'M';
S.montage = '/Users/christinadelta/Desktop/os_data/beads/spmdir/jobs/SPMeeg_montage.mat';
S.keepothers = 0;
S.keepsensors = 1;
S.updatehistory = 1;
D = spm_eeg_montage(S);


