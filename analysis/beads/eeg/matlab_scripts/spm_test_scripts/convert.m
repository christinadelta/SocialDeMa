spm('defaults', 'eeg');

%% define dataset S
S = [];
S.dataset = '/Users/christinadelta/Desktop/os_data/beads/subs/sub-01/sub_01_beads_block_01.bdf';
S.mode = 'continuous';
%%
S.channels = {
              'Fp1'
              'AF7'
              'AF3'
              'F1'
              'F3'
              'F5'
              'F7'
              'FT7'
              'FC5'
              'FC3'
              'FC1'
              'C1'
              'C3'
              'C5'
              'T7'
              'TP7'
              'CP5'
              'CP3'
              'CP1'
              'P1'
              'P3'
              'P5'
              'P7'
              'P9'
              'PO7'
              'PO3'
              'O1'
              'Iz'
              'Oz'
              'POz'
              'Pz'
              'CPz'
              'Fpz'
              'Fp2'
              'AF8'
              'AF4'
              'AFz'
              'Fz'
              'F2'
              'F4'
              'F6'
              'F8'
              'FT8'
              'FC6'
              'FC4'
              'FC2'
              'FCz'
              'Cz'
              'C2'
              'C4'
              'C6'
              'T8'
              'TP8'
              'CP6'
              'CP4'
              'CP2'
              'P2'
              'P4'
              'P6'
              'P8'
              'P10'
              'PO8'
              'PO4'
              'O2'
              'EXG3'
              'EXG4'
              'EXG5'
              'EXG6'
              }';
%% convert dataset S to spm object D
S.eventpadding = 0;
S.blocksize = 3276800;
S.checkboundary = 1;
S.saveorigheader = 0;
S.outfile = 'spmeeg_sub_01_beads_block_01';
S.timewin = [];
S.conditionlabels = {'Undefined'};
S.inputformat = [];
D = spm_eeg_convert(S);


