spm('defaults', 'eeg');

S = [];
S.D = '/Users/christinadelta/Desktop/os_data/beads/spmdir/output_tfr/rmtf_cefdfMspmeeg_sub_01_beads_block_01.mat';
S.c = [-1 1 0 0];
S.label = {'diff_difference'};
S.weighted = 1;
S.prefix = 'w';
D = spm_eeg_contrast(S);


