spm('defaults', 'eeg');

S = [];
S.D = '/Users/christinadelta/Desktop/os_data/beads/spmdir/output/spmeeg_sub_01_beads_block_01.mat';
S.mode = 'write';
S.blocksize = 655360;
S.prefix = 'M';
S.montage = '/Users/christinadelta/Desktop/os_data/beads/spmdir/jobs/SPMeeg_montage_block1.mat';
S.keepothers = 0;
S.keepsensors = 1;
S.updatehistory = 1;
D = spm_eeg_montage(S);


S = [];
S.D = '/Users/christinadelta/Desktop/os_data/beads/spmdir/output_erps/efdfMspmeeg_sub_01_beads_block_01.mat';
S.mode = 'reject';
S.badchanthresh = 0.2;
S.prefix = 'a';
S.append = true;
S.methods.channels = {'EEG'};
S.methods.fun = 'threshchan';
S.methods.settings.threshold = 100;
S.methods.settings.excwin = 1000;
D = spm_eeg_artefact(S);


