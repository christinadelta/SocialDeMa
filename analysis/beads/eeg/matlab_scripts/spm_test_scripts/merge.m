spm('defaults', 'eeg');

S = [];
S.D = [
       '/Users/christinadelta/Desktop/os_data/beads/spmdir/output/befdfMspmeeg_sub_01_beads_block_01.mat'
       '/Users/christinadelta/Desktop/os_data/beads/spmdir/output/befdfMspmeeg_sub_01_beads_block_02.mat'
       '/Users/christinadelta/Desktop/os_data/beads/spmdir/output/befdfMspmeeg_sub_01_beads_block_03.mat'
       '/Users/christinadelta/Desktop/os_data/beads/spmdir/output/befdfMspmeeg_sub_01_beads_block_04.mat'
       ];
S.recode.file = '.*';
S.recode.labelorg = '.*';
S.recode.labelnew = '#labelorg#';
S.prefix = 'c';
D = spm_eeg_merge(S);


