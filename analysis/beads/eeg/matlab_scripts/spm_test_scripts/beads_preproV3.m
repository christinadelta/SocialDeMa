%% PREPROCESSING AND ERP ANALYSIS

%% init spm 

clear all
clc

spm('defaults', 'eeg');

%% create S struct for conversion 
S = [];
S.dataset = '/Users/christinadelta/Desktop/os_data/beads/subs/sub-01/sub_01_beads_block_01.bdf';
% S.dataset = '/Users/christinadelta/Desktop/os_data/beads/subs/sub-01/sub_01_beads_block_02.bdf';
% S.dataset = '/Users/christinadelta/Desktop/os_data/beads/subs/sub-01/sub_01_beads_block_03.bdf';
% S.dataset = '/Users/christinadelta/Desktop/os_data/beads/subs/sub-01/sub_01_beads_block_04.bdf';
S.mode = 'continuous';

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
              'EXG1'
              'EXG2'
              'EXG3'
              'EXG4'
              'EXG5'
              'EXG6'
              }';

% convert bdf file to spm object and D struct
S.eventpadding = 0;
S.blocksize = 3276800;
S.checkboundary = 1;
S.saveorigheader = 0;
S.outfile = 'spmeeg_sub_01_beads_block_01';
% S.outfile = 'spmeeg_sub_01_beads_block_02';
% S.outfile = 'spmeeg_sub_01_beads_block_03';
% S.outfile = 'spmeeg_sub_01_beads_block_04';
S.timewin = [];
S.conditionlabels = {'Undefined'};
S.inputformat = [];
D = spm_eeg_convert(S);

%% create montage 
S = [];
% S.D = '/Users/christinadelta/Desktop/os_data/beads/spmdir/output_erps/spmeeg_sub_01_beads_block_01.mat';
% S.D = '/Users/christinadelta/Desktop/os_data/beads/spmdir/output_erps/spmeeg_sub_01_beads_block_02.mat';
% S.D = '/Users/christinadelta/Desktop/os_data/beads/spmdir/output_erps/spmeeg_sub_01_beads_block_03.mat';
S.D = '/Users/christinadelta/Desktop/os_data/beads/spmdir/output_erps/spmeeg_sub_01_beads_block_04.mat';

% Run Montage creation function beads_montage
% This function re-references by averaging across all electrodes. However,
% an initial step requires knowing in advance if there is a noisy channel
% and exclude it from averaging. Montage creation thus, need to be done
% individually for every subject (e.g. sub 1 has one noisy channel [channel
% 25]. This needs to be removed from averaging when re-referencing 

S.mode = 'write';
S.blocksize = 655360;
S.prefix = 'M';
S.montage = '/Users/christinadelta/Desktop/os_data/beads/spmdir/jobs/SPMeeg_montage_block4.mat';
S.keepothers = 0;
S.keepsensors = 1;
S.updatehistory = 1;
D = spm_eeg_montage(S);

%% high-pass filter 
S = [];
% S.D = '/Users/christinadelta/Desktop/os_data/beads/spmdir/output_erps/Mspmeeg_sub_01_beads_block_01.mat';
% S.D = '/Users/christinadelta/Desktop/os_data/beads/spmdir/output_erps/Mspmeeg_sub_01_beads_block_02.mat';
% S.D = '/Users/christinadelta/Desktop/os_data/beads/spmdir/output_erps/Mspmeeg_sub_01_beads_block_03.mat';
S.D = '/Users/christinadelta/Desktop/os_data/beads/spmdir/output_erps/Mspmeeg_sub_01_beads_block_04.mat';
S.type = 'butterworth';
S.band = 'high';
S.freq = 0.5;
S.dir = 'twopass';
S.order = 5;
S.prefix = 'f';
D = spm_eeg_filter(S);

%% downsample
S = [];
% S.D = '/Users/christinadelta/Desktop/os_data/beads/spmdir/output_erps/fMspmeeg_sub_01_beads_block_01.mat';
% S.D = '/Users/christinadelta/Desktop/os_data/beads/spmdir/output_erps/fMspmeeg_sub_01_beads_block_02.mat';
% S.D = '/Users/christinadelta/Desktop/os_data/beads/spmdir/output_erps/fMspmeeg_sub_01_beads_block_03.mat';
S.D = '/Users/christinadelta/Desktop/os_data/beads/spmdir/output_erps/fMspmeeg_sub_01_beads_block_04.mat';
S.fsample_new = 256;
S.method = 'resample';
S.prefix = 'd';
D = spm_eeg_downsample(S);

%% low-pass filter 
S = [];
% S.D = '/Users/christinadelta/Desktop/os_data/beads/spmdir/output_erps/dfMspmeeg_sub_01_beads_block_01.mat';
% S.D = '/Users/christinadelta/Desktop/os_data/beads/spmdir/output_erps/dfMspmeeg_sub_01_beads_block_02.mat';
% S.D = '/Users/christinadelta/Desktop/os_data/beads/spmdir/output_erps/dfMspmeeg_sub_01_beads_block_03.mat';
S.D = '/Users/christinadelta/Desktop/os_data/beads/spmdir/output_erps/dfMspmeeg_sub_01_beads_block_04.mat';
S.type = 'butterworth';
S.band = 'low';
S.freq = 30;
S.dir = 'twopass';
S.order = 5;
S.prefix = 'f';

D = spm_eeg_filter(S);

%% epoch data 

% FIRST SPECIFY TRIALS -- CREATE TRIAL DEFINITION MAT FILE
% this part creates a mat file (trial definition) using the the previously
% saved structure (low-pass filtered). This mat file should contain:
% 1. source of data
% 2. time window
% 3. trialdef (condition labels, event types, event values, trl shift)
% 4. condition labels (cell with conditions/strings)
% 5. trl (trial start, trial end, offset)

% use trial definition mat file for the S structure (to be converted into
% epoched MEEG object)
S = [];
% S.D = '/Users/christinadelta/Desktop/os_data/beads/spmdir/output_erps/fdfMspmeeg_sub_01_beads_block_01.mat';
% S.D = '/Users/christinadelta/Desktop/os_data/beads/spmdir/output_erps/fdfMspmeeg_sub_01_beads_block_02.mat';
% S.D = '/Users/christinadelta/Desktop/os_data/beads/spmdir/output_erps/fdfMspmeeg_sub_01_beads_block_03.mat';
S.D = '/Users/christinadelta/Desktop/os_data/beads/spmdir/output_erps/fdfMspmeeg_sub_01_beads_block_04.mat';
S.trialdef(1).conditionlabel = 'easydraw';
S.trialdef(1).eventtype = 'STATUS';
S.trialdef(1).eventvalue = 1;
S.trialdef(1).trlshift = 0;
S.trialdef(2).conditionlabel = 'easyurn';
S.trialdef(2).eventtype = 'STATUS';
S.trialdef(2).eventvalue = 2;
S.trialdef(2).trlshift = 0;
S.trialdef(3).conditionlabel = 'diffdraw';
S.trialdef(3).eventtype = 'STATUS';
S.trialdef(3).eventvalue = 3;
S.trialdef(3).trlshift = 0;
S.trialdef(4).conditionlabel = 'diffurn';
S.trialdef(4).eventtype = 'STATUS';
S.trialdef(4).eventvalue = 4;
S.trialdef(4).trlshift = 0;
S.timewin = [-500
             800];

% run the trial definition function here to get the trl matrix and
% condition labels
[trl, conditionlabels, S] = beads_trialdef(S);

S.bc = 1;
S.prefix = 'e';
S.eventpadding = 0;
D = spm_eeg_epochs(S);

% %% baseline correct 
% S = [];
% S.D = '/Users/christinadelta/Desktop/os_data/beads/spmdir/output/efdfMspmeeg_sub_01_beads_block_01.mat';
% % S.D = '/Users/christinadelta/Desktop/os_data/beads/spmdir/output/efdfMspmeeg_sub_01_beads_block_02.mat';
% %S.D = '/Users/christinadelta/Desktop/os_data/beads/spmdir/output/efdfMspmeeg_sub_01_beads_block_03.mat';
% % S.D = '/Users/christinadelta/Desktop/os_data/beads/spmdir/output/efdfMspmeeg_sub_01_beads_block_04.mat';
% S.timewin = [-500 -50];
% S.prefix = 'b';
% S.save = 1;
% S.updatehistory = 1;
% D = spm_eeg_bc(S);

%% preprocessing of block data should end here. 

% the next step is merging the blocks 

%% merge all blocks/runs 

S = [];
S.D = [
       '/Users/christinadelta/Desktop/os_data/beads/spmdir/output_erps/efdfMspmeeg_sub_01_beads_block_01.mat'
       '/Users/christinadelta/Desktop/os_data/beads/spmdir/output_erps/efdfMspmeeg_sub_01_beads_block_02.mat'
       '/Users/christinadelta/Desktop/os_data/beads/spmdir/output_erps/efdfMspmeeg_sub_01_beads_block_03.mat'
       '/Users/christinadelta/Desktop/os_data/beads/spmdir/output_erps/efdfMspmeeg_sub_01_beads_block_04.mat'
       ];
S.recode.file = '.*';
S.recode.labelorg = '.*';
S.recode.labelnew = '#labelorg#';
S.prefix = 'c';
D = spm_eeg_merge(S);


%% define cordinates/channel locations 
S = [];
S.D = '/Users/christinadelta/Desktop/os_data/beads/spmdir/output_erps/cefdfMspmeeg_sub_01_beads_block_01.mat';
S.task = 'defaulteegsens';
S.save = 1;
D = spm_eeg_prep(S);


S = [];
S.task = 'setcoor2d';
S.D = '/Users/christinadelta/Desktop/os_data/beads/spmdir/output_erps/cefdfMspmeeg_sub_01_beads_block_01.mat';
S.xy = [0.358127431070432 0.230361291206221 0.364687127498288 0.412094935710624 0.319499558310764 0.220586712887751 0.121753409066311 0.0551811363296929 0.178452330461275 0.296722753529597 0.402418214252759 0.402232139889891 0.295640676509763 0.176774133670791 0.05 0.0952050213281216 0.204129030986282 0.31059947251063 0.407784984615512 0.419107072168501 0.337261391153001 0.255936320417754 0.176983444035271 0.0990531376307864 0.274933610032685 0.374805341575767 0.381476047145341 0.485128199634413 0.490099574981175 0.493654015274171 0.496027125746896 0.497486874012439 0.494127551616874 0.63115333565553 0.764761267196509 0.633966624967236 0.497096812715531 0.498387742743428 0.590119958154716 0.684446454131226 0.784927979314371 0.878951869887076 0.947194362550941 0.826079267890216 0.706962412478808 0.597352645199803 0.498716216854912 0.498366895939448 0.598675341202977 0.70801035314695 0.826759976212305 0.95 0.898334363412286 0.795911344188124 0.69078343536388 0.593172800977718 0.57877084537302 0.658765804468392 0.735389719899426 0.810437806844084 0.879007062596769 0.708438648401535 0.611863922610018 0.599476237529535
        0.94291962707719 0.889006557470294 0.844625065001036 0.711929958683174 0.727189162575643 0.752686038230401 0.785972668551047 0.638824390760336 0.622237843107835 0.612055755800727 0.604719623778926 0.504717514961853 0.499865625106497 0.489925660632489 0.477921883184179 0.346661648796377 0.375837233637139 0.399737181066574 0.41403150193654 0.321035174889695 0.30567474167612 0.277954382109169 0.241048698933626 0.191392522403548 0.173810943732217 0.218490107582391 0.139058031510229 0.05 0.134825637229592 0.229144765608567 0.325406282673219 0.415924315197709 0.95 0.947839059007246 0.895554940224051 0.844643723365574 0.824931733550042 0.708258571334129 0.714599145507353 0.73251653002839 0.754709519350322 0.7891126332849 0.637244898999659 0.620328254919284 0.610930378565758 0.604973846888264 0.602533725218943 0.506171957423089 0.501865455408111 0.495100231864275 0.484556343615868 0.468273841263323 0.33147104080173 0.369228471090915 0.394026914867007 0.410294741704492 0.320570318800589 0.301276182177249 0.266698929482614 0.226857688879681 0.170533558640877 0.16331740311369 0.210477257752421 0.133428670630485];
%%
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

%save the file 
S.save = 1;
D = spm_eeg_prep(S);

%% artefact rejection

S = [];
S.D = '/Users/christinadelta/Desktop/os_data/beads/spmdir/output_erps/cefdfMspmeeg_sub_01_beads_block_01.mat';
S.mode = 'reject';
S.badchanthresh = 0.2;
S.prefix = 'a';
S.append = true;
S.methods.channels = {'EEG'};
S.methods.fun = 'threshchan';
S.methods.settings.threshold = 100;
S.methods.settings.excwin = 1000;
D = spm_eeg_artefact(S);

%% average conditions (Averaged ERPs)

S = [];
S.D = '/Users/christinadelta/Desktop/os_data/beads/spmdir/output_erps/acefdfMspmeeg_sub_01_beads_block_01.mat';
S.robust.ks = 3;
S.robust.bycondition = true;
S.robust.savew = false;
S.robust.removebad = false;
S.circularise = false;
S.prefix = 'm';
D = spm_eeg_average(S);


%% convert to 3D volums 

S = [];
% S.D = '/Users/christinadelta/Desktop/os_data/beads/spmdir/output_erps/macefdfMspmeeg_sub_01_beads_block_01.mat';
S.D = '/Users/christinadelta/Desktop/os_data/beads/spmdir/output_erps/acefdfMspmeeg_sub_01_beads_block_01.mat';
S.mode = 'scalp x time';
S.conditions = {};
S.channels = 'EEG';
S.timewin = [-Inf Inf];
S.freqwin = [-Inf Inf];
S.prefix = '';

D = spm_eeg_convert2images(S);

%% contrast ERP averaged conditions

%%% I will need to work on this %%%

matlabbatch{1}.spm.util.imcalc.input = {
                                        '/Users/christinadelta/Desktop/os_data/beads/spmdir/output/macbefdfMspmeeg_sub_01_beads_block_01/condition_easyurn.nii,1'
                                        '/Users/christinadelta/Desktop/os_data/beads/spmdir/output/macbefdfMspmeeg_sub_01_beads_block_01/condition_easydraw.nii,1'
                                        };
matlabbatch{1}.spm.util.imcalc.output = 'easy_output';
matlabbatch{1}.spm.util.imcalc.outdir = {''};
matlabbatch{1}.spm.util.imcalc.expression = 'i1 -i2';
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;

%% --------------------------------------------------------------------%%
% RUN TIME-FREQUENCY ANALYSIS
% write comments about tfr here

%% Time-Frequency Morlet Decomposition
S = [];
S.D = '/Users/christinadelta/Desktop/os_data/beads/spmdir/output_erps/cefdfMspmeeg_sub_01_beads_block_01.mat';
S.channels = {'all'};
S.frequencies = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55];
S.timewin = [-Inf Inf];
S.phase = 1;
S.method = 'morlet';
S.settings.ncycles = 7;
S.settings.timeres = 0;
S.settings.subsample = 5;
S.prefix = '';
D = spm_eeg_tf(S);

%% average power 
S = [];
S.D = '/Users/christinadelta/Desktop/os_data/beads/spmdir/output_tfr/tf_cefdfMspmeeg_sub_01_beads_block_01.mat';
S.robust = false;
S.circularise = false;
S.prefix = 'm';
D = spm_eeg_average(S);

%% average phase 
S = [];
S.D = '/Users/christinadelta/Desktop/os_data/beads/spmdir/output_tfr/tph_cefdfMspmeeg_sub_01_beads_block_01.mat';
S.robust = false;
S.circularise = true;
S.prefix = 'm';
D = spm_eeg_average(S);

%% baseline rescaling (only power file)

S = [];
S.D = '/Users/christinadelta/Desktop/os_data/beads/spmdir/output_tfr/mtf_cefdfMspmeeg_sub_01_beads_block_01.mat';
S.method = 'LogR';
S.prefix = 'r';
S.timewin = [-500 -50];
S.pooledbaseline = 0;
D = spm_eeg_tf_rescale(S);

%% contrast power and phase-locking files
S = [];
S.D = '/Users/christinadelta/Desktop/os_data/beads/spmdir/output_tfr/rmtf_cefdfMspmeeg_sub_01_beads_block_01.mat';
% S.D = '/Users/christinadelta/Desktop/os_data/beads/spmdir/output_tfr/mtph_cefdfMspmeeg_sub_01_beads_block_01.mat';
S.c = [-1 1 0 0]; % difficult difference 
% S.c = [0 0 -1 1]; % easy difference 
S.label = {'diff_difference'};
S.weighted = 1;
S.prefix = 'w';
D = spm_eeg_contrast(S);

%% convert TFR to images


%%  extract evoked (ERP) parietal data for analysis with behaviour and model

S = [];
S.D = '/Users/christinadelta/Desktop/os_data/beads/spmdir/output_erps/acefdfMspmeeg_sub_01_beads_block_01.mat';
S.timewin = [200 600];
S.freqwin = [-Inf Inf];

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

S.prefix = 'p';
D = spm_eeg_crop(S);

%% extract tfr parietal data (power and phase)

S = [];
% S.D =
% '/Users/christinadelta/Desktop/os_data/beads/spmdir/output_tfr/tf_cefdfMspmeeg_sub_01_beads_block_01.mat';
S.D = '/Users/christinadelta/Desktop/os_data/beads/spmdir/output_tfr/tph_cefdfMspmeeg_sub_01_beads_block_01.mat';
S.timewin = [200 600];
S.freqwin = [13 55];

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

S.prefix = 'p';
D = spm_eeg_crop(S);