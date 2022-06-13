spm('defaults', 'eeg');

S = [];
S.D = '/Users/christinadelta/Desktop/os_data/beads/spmdir/output_erps/cefdfMspmeeg_sub_01_beads_block_01.mat';
S.channels = {'all'};
S.frequencies = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50];
S.timewin = [-Inf Inf];
S.phase = 1;
S.method = 'morlet';
S.settings.ncycles = 7;
S.settings.timeres = 0;
S.settings.subsample = 5;
S.prefix = '';
D = spm_eeg_tf(S);


S = [];
S.D = '/Users/christinadelta/Desktop/os_data/beads/spmdir/output_tfr/tf_cefdfMspmeeg_sub_01_beads_block_01.mat';
S.robust = false;
S.circularise = false;
S.prefix = 'm';
D = spm_eeg_average(S);


S = [];
S.D = '/Users/christinadelta/Desktop/os_data/beads/spmdir/output_tfr/mtf_cefdfMspmeeg_sub_01_beads_block_01.mat';
S.method = 'LogR';
S.prefix = 'r';
S.timewin = [-500 -50];
S.pooledbaseline = 0;
D = spm_eeg_tf_rescale(S);


