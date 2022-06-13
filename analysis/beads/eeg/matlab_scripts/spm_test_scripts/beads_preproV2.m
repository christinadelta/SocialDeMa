spm('defaults', 'eeg');

S = [];
S.dataset = '/Users/christinadelta/Desktop/os_data/beads/subs/sub-01/sub_01_beads_block_01.bdf';
S.mode = 'continuous';
%% convert 
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

S.eventpadding = 0;
S.blocksize = 3276800;
S.checkboundary = 1;
S.saveorigheader = 0;
S.outfile = 'spmeeg_sub_01_beads_block_01';
S.timewin = [];
S.conditionlabels = {'Undefined'};
S.inputformat = [];
D = spm_eeg_convert(S);

%% montage
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

%% highpass filter
S = [];
S.D = '/Users/christinadelta/Desktop/os_data/beads/spmdir/output/Mspmeeg_sub_01_beads_block_01.mat';
S.type = 'butterworth';
S.band = 'high';
S.freq = 0.1;
S.dir = 'twopass';
S.order = 5;
S.prefix = 'f';
D = spm_eeg_filter(S);

%% resample
S = [];
S.D = '/Users/christinadelta/Desktop/os_data/beads/spmdir/output/fMspmeeg_sub_01_beads_block_01.mat';
S.fsample_new = 256;
S.method = 'resample';
S.prefix = 'd';
D = spm_eeg_downsample(S);

%% lowpass filter 
S = [];
S.D = '/Users/christinadelta/Desktop/os_data/beads/spmdir/output/dfMspmeeg_sub_01_beads_block_01.mat';
S.type = 'butterworth';
S.band = 'low';
S.freq = 30;
S.dir = 'twopass';
S.order = 5;
S.prefix = 'f';
D = spm_eeg_filter(S);

%% define trials/epochs
S = [];
S.D = '/Users/christinadelta/Desktop/os_data/beads/spmdir/output/fdfMspmeeg_sub_01_beads_block_01.mat';
S.trialdef(1).conditionlabel = 'start';
S.trialdef(1).eventtype = 'STATUS';
S.trialdef(1).eventvalue = 102;
S.trialdef(1).trlshift = 0;
S.trialdef(2).conditionlabel = 'end';
S.trialdef(2).eventtype = 'STATUS';
S.trialdef(2).eventvalue = 103;
S.trialdef(2).trlshift = 0;
S.trialdef(3).conditionlabel = 'easy';
S.trialdef(3).eventtype = 'STATUS';
S.trialdef(3).eventvalue = 1;
S.trialdef(3).trlshift = 0;
S.trialdef(4).conditionlabel = 'easy';
S.trialdef(4).eventtype = 'STATUS';
S.trialdef(4).eventvalue = 2;
S.trialdef(4).trlshift = 0;
S.trialdef(5).conditionlabel = 'difficult';
S.trialdef(5).eventtype = 'STATUS';
S.trialdef(5).eventvalue = 3;
S.trialdef(5).trlshift = 0;
S.trialdef(6).conditionlabel = 'difficult';
S.trialdef(6).eventtype = 'STATUS';
S.trialdef(6).eventvalue = 4;
S.trialdef(6).trlshift = 0;
S.timewin = [-500
             800];
S.bc = 0;
S.prefix = 'e';
S.eventpadding = 0;
D = spm_eeg_epochs(S);

%% baseline correction
S = [];
S.D = '/Users/christinadelta/Desktop/os_data/beads/spmdir/output/efdfMspmeeg_sub_01_beads_block_01.mat';
S.timewin = [-500 -50];
S.prefix = 'b';
S.save = 1;
S.updatehistory = 1;
D = spm_eeg_bc(S);

%% define coordinates/ channel connection
S = [];
S.D = '/Users/christinadelta/Desktop/os_data/beads/spmdir/output/befdfMspmeeg_sub_01_beads_block_01.mat';
S.task = 'defaulteegsens';
S.save = 1;
D = spm_eeg_prep(S);


S = [];
S.task = 'setcoor2d';
S.D = '/Users/christinadelta/Desktop/os_data/beads/spmdir/output/befdfMspmeeg_sub_01_beads_block_01.mat';
S.xy = [0.358127431070432 0.230361291206221 0.364687127498288 0.412094935710624 0.319499558310764 0.220586712887751 0.121753409066311 0.0551811363296929 0.178452330461275 0.296722753529597 0.402418214252759 0.402232139889891 0.295640676509763 0.176774133670791 0.05 0.0952050213281216 0.204129030986282 0.31059947251063 0.407784984615512 0.419107072168501 0.337261391153001 0.255936320417754 0.176983444035271 0.0990531376307864 0.274933610032685 0.374805341575767 0.381476047145341 0.485128199634413 0.490099574981175 0.493654015274171 0.496027125746896 0.497486874012439 0.494127551616874 0.63115333565553 0.764761267196509 0.633966624967236 0.497096812715531 0.498387742743428 0.590119958154716 0.684446454131226 0.784927979314371 0.878951869887076 0.947194362550941 0.826079267890216 0.706962412478808 0.597352645199803 0.498716216854912 0.498366895939448 0.598675341202977 0.70801035314695 0.826759976212305 0.95 0.898334363412286 0.795911344188124 0.69078343536388 0.593172800977718 0.57877084537302 0.658765804468392 0.735389719899426 0.810437806844084 0.879007062596769 0.708438648401535 0.611863922610018 0.599476237529535
        0.94291962707719 0.889006557470294 0.844625065001036 0.711929958683174 0.727189162575643 0.752686038230401 0.785972668551047 0.638824390760336 0.622237843107835 0.612055755800727 0.604719623778926 0.504717514961853 0.499865625106497 0.489925660632489 0.477921883184179 0.346661648796377 0.375837233637139 0.399737181066574 0.41403150193654 0.321035174889695 0.30567474167612 0.277954382109169 0.241048698933626 0.191392522403548 0.173810943732217 0.218490107582391 0.139058031510229 0.05 0.134825637229592 0.229144765608567 0.325406282673219 0.415924315197709 0.95 0.947839059007246 0.895554940224051 0.844643723365574 0.824931733550042 0.708258571334129 0.714599145507353 0.73251653002839 0.754709519350322 0.7891126332849 0.637244898999659 0.620328254919284 0.610930378565758 0.604973846888264 0.602533725218943 0.506171957423089 0.501865455408111 0.495100231864275 0.484556343615868 0.468273841263323 0.33147104080173 0.369228471090915 0.394026914867007 0.410294741704492 0.320570318800589 0.301276182177249 0.266698929482614 0.226857688879681 0.170533558640877 0.16331740311369 0.210477257752421 0.133428670630485];

    
S.label = {
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
           }';

S.save = 1;
D = spm_eeg_prep(S);

%% convert to images

matlabbatch{1}.spm.meeg.images.convert2images.D = {'/Users/christinadelta/Desktop/os_data/beads/spmdir/output/macbefdfMspmeeg_sub_01_beads_block_01.mat'};
matlabbatch{1}.spm.meeg.images.convert2images.mode = 'scalp x time';
matlabbatch{1}.spm.meeg.images.convert2images.conditions = {};
matlabbatch{1}.spm.meeg.images.convert2images.channels{1}.type = 'EEG';
matlabbatch{1}.spm.meeg.images.convert2images.timewin = [-Inf Inf];
matlabbatch{1}.spm.meeg.images.convert2images.freqwin = [-Inf Inf];
matlabbatch{1}.spm.meeg.images.convert2images.prefix = '';
