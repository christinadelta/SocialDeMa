spm('defaults', 'eeg');

S = [];
S.D = '/Users/christinadelta/Desktop/os_data/beads/spmdir/output_tfr/tf_cefdfMspmeeg_sub_01_beads_block_01.mat';
S.timewin = [200 600];
S.freqwin = [13 50];
%%
S.channels = {
              'P1'
              'P3'
              'P5'
              'P7'
              'P9'
              'Pz'
              'P2'
              'P4'
              'P6'
              'P8'
              'P10'
              }';
%%
S.prefix = 'p';
D = spm_eeg_crop(S);


