spm('defaults', 'eeg');

S = [];
S.D = '/Users/christinadelta/Desktop/os_data/beads/spmdir/output/abefdfMspmeeg_sub_01_beads_block_01.mat';
S.robust.ks = 3;
S.robust.bycondition = false;
S.robust.savew = false;
S.robust.removebad = false;
S.circularise = false;
S.prefix = 'm';
D = spm_eeg_average(S);


