spm('defaults', 'eeg');

S = [];
S.D = '/Users/christinadelta/Desktop/os_data/beads/spmdir/output/cbefdfMspmeeg_sub_01_beads_block_01.mat';
S.mode = 'reject';
S.badchanthresh = 0.2;
S.prefix = 'a';
S.append = true;
S.methods.channels = {'all'};
S.methods.fun = 'threshchan';
S.methods.settings.threshold = 500;
S.methods.settings.excwin = 1000;
D = spm_eeg_artefact(S);


S = [];
S.D = '/Users/christinadelta/Desktop/os_data/beads/spmdir/output/acbefdfMspmeeg_sub_01_beads_block_01.mat';
S.robust.ks = 3;
S.robust.bycondition = true;
S.robust.savew = false;
S.robust.removebad = false;
S.circularise = false;
S.prefix = 'm';
D = spm_eeg_average(S);


