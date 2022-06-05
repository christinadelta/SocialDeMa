spm('defaults', 'eeg');

%% define trials/epochs 
S = [];
S.D = '/Users/christinadelta/Desktop/os_data/beads/spmdir/output/fdfMspmeeg_sub_01_beads_block_01.mat';
S.trialdef(1).conditionlabel = 'easy';
S.trialdef(1).eventtype = 'STATUS';
S.trialdef(1).eventvalue = 1;
S.trialdef(1).trlshift = 0;
S.trialdef(2).conditionlabel = 'easy';
S.trialdef(2).eventtype = 'STATUS';
S.trialdef(2).eventvalue = 2;
S.trialdef(2).trlshift = 0;
S.trialdef(3).conditionlabel = 'difficult';
S.trialdef(3).eventtype = 'STATUS';
S.trialdef(3).eventvalue = 3;
S.trialdef(3).trlshift = 0;
S.trialdef(4).conditionlabel = 'difficult';
S.trialdef(4).eventtype = 'STATUS';
S.trialdef(4).eventvalue = 4;
S.trialdef(4).trlshift = 0;
S.timewin = [-500
             800];
S.bc = 1;
S.prefix = 'e';
S.eventpadding = 0;
D = spm_eeg_epochs(S);


