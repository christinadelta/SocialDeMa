% Script to visualise Beads EEG data 
% Created September 05/09/2022

% The script plots Beads EEG timeseries data for single (specific) sensors
% and for averaged sensors for codntitions: draw-urn choices (diff and
% easy).

%% Create required directories and define paths %%

% clear workspace
clear all
clc

% add paths
basedir     = pwd;
avDir       = fullfile(basedir, 'averages');
addpath(genpath(fullfile(basedir, 'scripts')));
addpath(genpath(fullfile(basedir, 'utilities')));

% init spm 
spm('defaults', 'eeg');

%% 1. PARIETAL - CROP AVERAGED EEG DATA FOR VISUALISATION %%

% 1. crop averaged ERP data for parietal electrodes only
% init S struct 
S               = [];
S.D             = fullfile(avDir, 'grand_average.mat');
S.timewin       = [-500 800];
S.freqwin       = [-Inf Inf];
S.channels      = {'P1','P3','P5','Pz','P2','P4','P6','P8'};
S.prefix        = 'par_';
D               = spm_eeg_crop(S);

% get access and extract the actual data
S               = [];
S.P             = fullfile(avDir, 'par_grand_average.mat');
D               = spm_eeg_load(S.P);

% access the data
data            = D(:,:,:);      % D(channels, samples, trials)
sub_conds       = conditions(D); % get conditions order for this sub

% load spm meeg with fieldtrip. This is to get acces to the samples (time points)
par_data        = spm2fieldtrip(D);
pardata_samples = par_data.time;

clear par_data


%% PLOT PARIETAL SENSORS

%%%% PLOT ALL PARIETAL SENSORS for diff and easy
% First seperate trials/conditions
diffdraw_par    = data(:,:,1); % 1. Difficult-draw trial
diffurn_par     = data(:,:,2); % 2. Difficult-urn trial
easydraw_par    = data(:,:,3); % 3. Easy-draw trial
easyurn_par     = data(:,:,4); % 4. Easy-urn trial

% average channels
diffdraw_av     = mean(diffdraw_par,1);
diffurn_av      = mean(diffurn_par,1);
easydraw_av     = mean(easydraw_par,1);
easyurn_av      = mean(easyurn_par,1);

%%% plot averaged parietal ERP for diff
plot(pardata_samples{1,2},diffdraw_av, '-b')
grid on % add labels 
xlabel('Time [s]')
ylabel('Field intensity [in uV]')
title('Averaged parietal sensors for 0.6 probability condition')
hold on % to also add cond b
plot(pardata_samples{1,2}, diffurn_av,'-r','MarkerSize',6),
legend('0.6-draw','0.6-urn','Location','northwest') % Add a legend in the upper left:

%%% plot averaged parietal ERP for easy
plot(pardata_samples{1,2},easydraw_av, '-b')
grid on % add labels 
xlabel('Time [s]')
ylabel('Field intensity [in uV]')
title('Averaged parietal sensors for 0.8 probability condition')
hold on % to also add cond b
plot(pardata_samples{1,2}, easyurn_av,'-r','MarkerSize',6),
legend('0.8-draw  ','0.8-urn  ','Location','northwest') % Add a legend in the upper left:

%%%% PLOT LEFT & RIGHT PARIETAL SENSORS for diff condition
% First seperate trials/conditions
diffdraw_parleft    = data(1:4,:,1); % 1. Difficult-draw left trial
diffdraw_paright    = data(5:8,:,1); % 2. Difficult-urn right trial
diffurn_parleft     = data(1:4,:,2); % 3. Difficult-urn left trial
diffurn_paright     = data(5:8,:,2); % 4. Difficult-urn right trial

% average over sensors
diffdraw_parleft_av = mean(diffdraw_parleft,1);
diffurn_parleft_av  = mean(diffurn_parleft,1);
diffdraw_paright_av = mean(diffdraw_paright,1);
diffurn_paright_av  = mean(diffurn_paright,1);

%%% plot averaged parietal ERP for 0.6 left sensors
plot(pardata_samples{1,2},diffdraw_parleft_av, '-b')
grid on % add labels 
xlabel('Time [s]')
ylabel('Field intensity [in uV]')
title('Averaged parietal left sensors for 0.6 probability condition')
hold on % to also add cond b
plot(pardata_samples{1,2}, diffurn_parleft_av,'-r','MarkerSize',6),
legend('0.6-draw  ','0.6-urn  ','Location','northwest') % Add a legend in the upper left:

%%% plot averaged parietal ERP for 0.6 right sensors
plot(pardata_samples{1,2},diffdraw_paright_av, '-b')
grid on % add labels 
xlabel('Time [s]')
ylabel('Field intensity [in uV]')
title('Averaged parietal right sensors for 0.6 probability condition')
hold on % to also add cond b
plot(pardata_samples{1,2}, diffurn_paright_av,'-r','MarkerSize',6),
legend('0.6-draw  ','0.6-urn  ','Location','northwest') % Add a legend in the upper left:

%%%% PLOT LEFT & RIGHT PARIETAL SENSORS for 0.8 condition
% First seperate trials/conditions
easydraw_parleft    = data(1:4,:,3); % 1. easy-draw left trial
easydraw_paright    = data(5:8,:,3); % 2. easy-draw right trial
easyurn_parleft     = data(1:4,:,4); % 3. wasy-urn left trial
easyurn_paright     = data(5:8,:,4); % 4. easy-urn right trial

% average over sensors
easydraw_parleft_av = mean(easydraw_parleft,1);
easyurn_parleft_av  = mean(easyurn_parleft,1);
easydraw_paright_av = mean(easydraw_paright,1);
easyurn_paright_av  = mean(easyurn_paright,1);

%%% plot averaged parietal ERP for 0.8 left sensors
plot(pardata_samples{1,2},easydraw_parleft_av, '-b')
grid on % add labels 
xlabel('Time [s]')
ylabel('Field intensity [in uV]')
title('Averaged parietal left sensors for 0.8 probability condition')
hold on % to also add cond b
plot(pardata_samples{1,2}, easyurn_parleft_av,'-r','MarkerSize',6),
legend('0.8-draw  ','0.8-urn  ','Location','northwest') % Add a legend in the upper left:

%%% plot averaged parietal ERP for 0.8 right sensors
plot(pardata_samples{1,2},easydraw_paright_av, '-b')
grid on % add labels 
xlabel('Time [s]')
ylabel('Field intensity [in uV]')
title('Averaged parietal right sensors for 0.8 probability condition')
hold on % to also add cond b
plot(pardata_samples{1,2}, easyurn_paright_av,'-r','MarkerSize',6),
legend('0.8-draw  ','0.8-urn  ','Location','northwest') % Add a legend in the upper left:

% clear stuff
clear diffdraw_av diffdraw_par diffdraw_paright diffdraw_paright_av diffdraw_parleft diffdraw_parleft_av data diffurn_av 
clear diffurn_par diffurn_paright diffurn_paright_av diffurn_parleft diffurn_parleft_av pardata_samples
clear easydraw_av easydraw_par easydraw_paright easydraw_paright_av easydraw_parleft easydraw_parleft_av
clear easyurn_av easyurn_par easyurn_paright easyurn_paright_av easyurn_parleft_av easyurn_parleft sub_conds D S

%% 2. FRONTAL - CROP AVERAGED EEG DATA FOR VISUALISATION %%

% 1. crop averaged ERP data for parietal electrodes only
% init S struct 
S               = [];
S.D             = fullfile(avDir, 'grand_average.mat');
S.timewin       = [-500 800];
S.freqwin       = [-Inf Inf];
S.channels      = {'F1','F3','F5','F7','Fz','F2','F4','F6','F8'};
S.prefix        = 'front_';
D               = spm_eeg_crop(S);

% get access and extract the actual data
S               = [];
S.P             = fullfile(avDir, 'front_grand_average.mat');
D               = spm_eeg_load(S.P);

% access the data
data            = D(:,:,:);      % D(channels, samples, trials)
sub_conds       = conditions(D); % get conditions order for this sub

% load spm meeg with fieldtrip. This is to get acces to the samples (time points)
front_data        = spm2fieldtrip(D);
frontdata_samples = front_data.time;

clear front_data

%% PLOT FRONTAL SENSORS

%%%% PLOT ALL FRONTAL SENSORS for diff and easy
% First seperate trials/conditions
diffdraw_fr    = data(:,:,1); % 1. Difficult-draw trial
diffurn_fr     = data(:,:,2); % 2. Difficult-urn trial
easydraw_fr    = data(:,:,3); % 3. Easy-draw trial
easyurn_fr     = data(:,:,4); % 4. Easy-urn trial

% average channels
diffdraw_av     = mean(diffdraw_fr,1);
diffurn_av      = mean(diffurn_fr,1);
easydraw_av     = mean(easydraw_fr,1);
easyurn_av      = mean(easyurn_fr,1);

%%% plot averaged frontal ERP for diff
plot(frontdata_samples{1,2},diffdraw_av, '-b')
grid on % add labels 
xlabel('Time [s]')
ylabel('Field intensity [in uV]')
title('Averaged frontal sensors for 0.6 probability condition')
hold on % to also add cond b
plot(frontdata_samples{1,2}, diffurn_av,'-r','MarkerSize',6),
legend('0.6-draw','0.6-urn','Location','northwest') % Add a legend in the upper left:

%%% plot averaged frontal ERP for easy
plot(frontdata_samples{1,2},easydraw_av, '-b')
grid on % add labels 
xlabel('Time [s]')
ylabel('Field intensity [in uV]')
title('Averaged frontal sensors for 0.8 probability condition')
hold on % to also add cond b
plot(frontdata_samples{1,2}, easyurn_av,'-r','MarkerSize',6),
legend('0.8-draw  ','0.8-urn  ','Location','northwest') % Add a legend in the upper left:

%%%% PLOT LEFT & RIGHT FRONTAL SENSORS for diff condition
% First seperate trials/conditions
diffdraw_frleft    = data(1:5,:,1); % 1. Difficult-draw left trial
diffdraw_fright    = data(6:9,:,1); % 2. Difficult-draw right trial
diffurn_frleft     = data(1:5,:,2); % 3. Difficult-urn left trial
diffurn_fright     = data(6:9,:,2); % 4. Difficult-urn right trial

% average over sensors
diffdraw_frleft_av = mean(diffdraw_frleft,1);
diffurn_frleft_av  = mean(diffurn_frleft,1);
diffdraw_fright_av = mean(diffdraw_fright,1);
diffurn_fright_av  = mean(diffurn_fright,1);

%%% plot averaged frontal ERP for 0.6 left sensors
plot(frontdata_samples{1,2},diffdraw_frleft_av, '-b')
grid on % add labels 
xlabel('Time [s]')
ylabel('Field intensity [in uV]')
title('Averaged frontal left sensors for 0.6 probability condition')
hold on % to also add cond b
plot(frontdata_samples{1,2}, diffurn_frleft_av,'-r','MarkerSize',6),
legend('0.6-draw  ','0.6-urn  ','Location','northwest') % Add a legend in the upper left:

%%% plot averaged frontal ERP for 0.6 right sensors
plot(frontdata_samples{1,2},diffdraw_fright_av, '-b')
grid on % add labels 
xlabel('Time [s]')
ylabel('Field intensity [in uV]')
title('Averaged frontal right sensors for 0.6 probability condition')
hold on % to also add cond b
plot(frontdata_samples{1,2}, diffurn_fright_av,'-r','MarkerSize',6),
legend('0.6-draw  ','0.6-urn  ','Location','northwest') % Add a legend in the upper left:

%%%% PLOT LEFT & RIGHT FRONTAL SENSORS for 0.8 condition
% First seperate trials/conditions
easydraw_frleft    = data(1:5,:,3); % 1. easy-draw left trial
easydraw_fright    = data(6:9,:,3); % 2. easy-draw right trial
easyurn_frleft     = data(1:5,:,4); % 3. wasy-urn left trial
easyurn_fright     = data(6:9,:,4); % 4. easy-urn right trial

% average over sensors
easydraw_frleft_av = mean(easydraw_frleft,1);
easyurn_frleft_av  = mean(easyurn_frleft,1);
easydraw_fright_av = mean(easydraw_fright,1);
easyurn_fright_av  = mean(easyurn_fright,1);

%%% plot averaged parietal ERP for 0.8 left sensors
plot(frontdata_samples{1,2},easydraw_frleft_av, '-b')
grid on % add labels 
xlabel('Time [s]')
ylabel('Field intensity [in uV]')
title('Averaged frontal left sensors for 0.8 probability condition')
hold on % to also add cond b
plot(frontdata_samples{1,2}, easyurn_frleft_av,'-r','MarkerSize',6),
legend('0.8-draw  ','0.8-urn  ','Location','northwest') % Add a legend in the upper left:

%%% plot averaged parietal ERP for 0.8 right sensors
plot(frontdata_samples{1,2},easydraw_fright_av, '-b')
grid on % add labels 
xlabel('Time [s]')
ylabel('Field intensity [in uV]')
title('Averaged frontal right sensors for 0.8 probability condition')
hold on % to also add cond b
plot(frontdata_samples{1,2}, easyurn_fright_av,'-r','MarkerSize',6),
legend('0.8-draw  ','0.8-urn  ','Location','northwest') % Add a legend in the upper left:

clear data diffdraw_av diffdraw_fr diffdraw_fright diffdraw_fright_av diffdraw_frleft diffdraw_frleft_av 
clear sub_conds S D frontdata_samples easydraw_av easydraw_fr easydraw_fright easydraw_fright_av 
clear easydraw_frleft easydraw_frleft_av easydraw_fr easydraw_fright easyurn_av easyurn_fr easyurn_fright_av diffurn_frleft
clear easyurn_fright easyurn_frleft easyurn_frleft_av diffurn_fr diffurn_av diffurn_fright diffurn_fright_av diffurn_frleft_av

%% 3. CENTRO-PARIETAL - CROP AVERAGED EEG DATA FOR VISUALISATION %%

% 1. crop averaged ERP data for parietal electrodes only
% init S struct 
S               = [];
S.D             = fullfile(avDir, 'grand_average.mat');
S.timewin       = [-500 800];
S.freqwin       = [-Inf Inf];
S.channels      = {'CP5','CP3','CP1','CPz','CP6','CP4','CP2'};
S.prefix        = 'cp_';
D               = spm_eeg_crop(S);

% get access and extract the actual data
S               = [];
S.P             = fullfile(avDir, 'cp_grand_average.mat');
D               = spm_eeg_load(S.P);

% access the data
data            = D(:,:,:);      % D(channels, samples, trials)
sub_conds       = conditions(D); % get conditions order for this sub

% load spm meeg with fieldtrip. This is to get acces to the samples (time points)
cp_data         = spm2fieldtrip(D);
cpdata_samples  = cp_data.time;

clear cp_data

%% PLOT PARIETAL SENSORS

%%%% PLOT ALL PARIETAL SENSORS for diff and easy
% First seperate trials/conditions
diffdraw_cp    = data(:,:,1); % 1. Difficult-draw trial
diffurn_cp     = data(:,:,2); % 2. Difficult-urn trial
easydraw_cp    = data(:,:,3); % 3. Easy-draw trial
easyurn_cp     = data(:,:,4); % 4. Easy-urn trial

% average channels
diffdraw_av     = mean(diffdraw_cp,1);
diffurn_av      = mean(diffurn_cp,1);
easydraw_av     = mean(easydraw_cp,1);
easyurn_av      = mean(easyurn_cp,1);

%%% plot averaged parietal ERP for diff
plot(cpdata_samples{1,2},diffdraw_av, '-b')
grid on % add labels 
xlabel('Time [s]')
ylabel('Field intensity [in uV]')
title('Averaged centro-parietal sensors for 0.6 probability condition')
hold on % to also add cond b
plot(cpdata_samples{1,2}, diffurn_av,'-r','MarkerSize',6),
legend('0.6-draw','0.6-urn','Location','northwest') % Add a legend in the upper left:

%%% plot averaged parietal ERP for easy
plot(cpdata_samples{1,2},easydraw_av, '-b')
grid on % add labels 
xlabel('Time [s]')
ylabel('Field intensity [in uV]')
title('Averaged centro-parietal sensors for 0.8 probability condition')
hold on % to also add cond b
plot(cpdata_samples{1,2}, easyurn_av,'-r','MarkerSize',6),
legend('0.8-draw  ','0.8-urn  ','Location','northwest') % Add a legend in the upper left:

%%%% PLOT LEFT & RIGHT PARIETAL SENSORS for diff condition
% First seperate trials/conditions
diffdraw_cpleft     = data(1:4,:,1); % 1. Difficult-draw left trial
diffdraw_cpright    = data(5:7,:,1); % 2. Difficult-draw right trial
diffurn_cpleft      = data(1:4,:,2); % 3. Difficult-urn left trial
diffurn_cpright     = data(5:7,:,2); % 4. Difficult-urn right trial

% average over sensors
diffdraw_cpleft_av = mean(diffdraw_cpleft,1);
diffurn_cpleft_av  = mean(diffurn_cpleft,1);
diffdraw_cpright_av = mean(diffdraw_cpright,1);
diffurn_cpright_av  = mean(diffurn_cpright,1);

%%% plot averaged parietal ERP for 0.6 left sensors
plot(cpdata_samples{1,2},diffdraw_cpleft_av, '-b')
grid on % add labels 
xlabel('Time [s]')
ylabel('Field intensity [in uV]')
title('Averaged centro-parietal left sensors for 0.6 probability condition')
hold on % to also add cond b
plot(cpdata_samples{1,2}, diffurn_cpleft_av,'-r','MarkerSize',6),
legend('0.6-draw  ','0.6-urn  ','Location','northwest') % Add a legend in the upper left:

%%% plot averaged parietal ERP for 0.6 right sensors
plot(cpdata_samples{1,2},diffdraw_cpright_av, '-b')
grid on % add labels 
xlabel('Time [s]')
ylabel('Field intensity [in uV]')
title('Averaged centro-parietal right sensors for 0.6 probability condition')
hold on % to also add cond b
plot(cpdata_samples{1,2}, diffurn_cpright_av,'-r','MarkerSize',6),
legend('0.6-draw  ','0.6-urn  ','Location','northwest') % Add a legend in the upper left:

%%%% PLOT LEFT & RIGHT PARIETAL SENSORS for 0.8 condition
% First seperate trials/conditions
easydraw_cpleft     = data(1:4,:,3); % 1. easy-draw left trial
easydraw_cpright    = data(5:7,:,3); % 2. easy-draw right trial
easyurn_cpleft      = data(1:4,:,4); % 3. wasy-urn left trial
easyurn_cpright     = data(5:7,:,4); % 4. easy-urn right trial

% average over sensors
easydraw_cpleft_av  = mean(easydraw_cpleft,1);
easyurn_cpleft_av   = mean(easyurn_cpleft,1);
easydraw_cpright_av = mean(easydraw_cpright,1);
easyurn_cpright_av  = mean(easyurn_cpright,1);

%%% plot averaged parietal ERP for 0.8 left sensors
plot(cpdata_samples{1,2},easydraw_cpleft_av, '-b')
grid on % add labels 
xlabel('Time [s]')
ylabel('Field intensity [in uV]')
title('Averaged centro-parietal left sensors for 0.8 probability condition')
hold on % to also add cond b
plot(cpdata_samples{1,2}, easyurn_cpleft_av,'-r','MarkerSize',6),
legend('0.8-draw  ','0.8-urn  ','Location','northwest') % Add a legend in the upper left:

%%% plot averaged parietal ERP for 0.8 right sensors
plot(cpdata_samples{1,2},easydraw_cpright_av, '-b')
grid on % add labels 
xlabel('Time [s]')
ylabel('Field intensity [in uV]')
title('Averaged centro-parietal right sensors for 0.8 probability condition')
hold on % to also add cond b
plot(cpdata_samples{1,2}, easyurn_cpright_av,'-r','MarkerSize',6),
legend('0.8-draw  ','0.8-urn  ','Location','northwest') % Add a legend in the upper left:

clear data D cpdata_samples diffdraw_av diffdraw_cp diffdraw_cpleft diffdraw_cpleft_av diffdraw_cpright diffdraw_cpright_av
clear diffurn_av diffurn_cp diffurn_cpleft diffurn_cpleft_av diffurn_cpright diffurn_cpright_av easydraw_av easydraw_cp
clear easydraw_cpleft easydraw_cpleft_av easyurn_av easyurn_cp easyurn_cpleft easyurn_cpleft_av easyurn_cpright 
clear easyurn_cpright_av S sub_conds easydraw_cpright_av easydraw_cpright

%% 4. FRONTO-CENTRAL - CROP AVERAGED EEG DATA FOR VISUALISATION %%

% 1. crop averaged ERP data for parietal electrodes only
% init S struct 
S               = [];
S.D             = fullfile(avDir, 'grand_average.mat');
S.timewin       = [-500 800];
S.freqwin       = [-Inf Inf];
S.channels      = {'FC5','FC3','FC1','FCz','FC6','FC4','FC2'};
S.prefix        = 'fc_';
D               = spm_eeg_crop(S);

% get access and extract the actual data
S               = [];
S.P             = fullfile(avDir, 'fc_grand_average.mat');
D               = spm_eeg_load(S.P);

% access the data
data            = D(:,:,:);      % D(channels, samples, trials)
sub_conds       = conditions(D); % get conditions order for this sub

% load spm meeg with fieldtrip. This is to get acces to the samples (time points)
fc_data         = spm2fieldtrip(D);
fcdata_samples  = fc_data.time;

clear fc_data

%% PLOT FRONTO-CENTRAL SENSORS

%%%% PLOT ALL PARIETAL SENSORS for diff and easy
% First seperate trials/conditions
diffdraw_fc    = data(:,:,1); % 1. Difficult-draw trial
diffurn_fc     = data(:,:,2); % 2. Difficult-urn trial
easydraw_fc    = data(:,:,3); % 3. Easy-draw trial
easyurn_fc     = data(:,:,4); % 4. Easy-urn trial

% average channels
diffdraw_av     = mean(diffdraw_fc,1);
diffurn_av      = mean(diffurn_fc,1);
easydraw_av     = mean(easydraw_fc,1);
easyurn_av      = mean(easyurn_fc,1);

%%% plot averaged parietal ERP for diff
plot(fcdata_samples{1,2},diffdraw_av, '-b')
grid on % add labels 
xlabel('Time [s]')
ylabel('Field intensity [in uV]')
title('Averaged fronto-central sensors for 0.6 probability condition')
hold on % to also add cond b
plot(fcdata_samples{1,2}, diffurn_av,'-r','MarkerSize',6),
legend('0.6-draw','0.6-urn','Location','northwest') % Add a legend in the upper left:

%%% plot averaged parietal ERP for easy
plot(fcdata_samples{1,2},easydraw_av, '-b')
grid on % add labels 
xlabel('Time [s]')
ylabel('Field intensity [in uV]')
title('Averaged fronto-central sensors for 0.8 probability condition')
hold on % to also add cond b
plot(fcdata_samples{1,2}, easyurn_av,'-r','MarkerSize',6),
legend('0.8-draw  ','0.8-urn  ','Location','northwest') % Add a legend in the upper left:

%%%% PLOT LEFT & RIGHT PARIETAL SENSORS for diff condition
% First seperate trials/conditions
diffdraw_fcleft     = data(1:4,:,1); % 1. Difficult-draw left trial
diffdraw_fcright    = data(5:7,:,1); % 2. Difficult-urn right trial
diffurn_fcleft      = data(1:4,:,2); % 3. Difficult-urn left trial
diffurn_fcright     = data(5:7,:,2); % 4. Difficult-urn right trial

% average over sensors
diffdraw_fcleft_av  = mean(diffdraw_fcleft,1);
diffurn_fcleft_av   = mean(diffurn_fcleft,1);
diffdraw_fcright_av = mean(diffdraw_fcright,1);
diffurn_fcright_av  = mean(diffurn_fcright,1);

%%% plot averaged parietal ERP for 0.6 left sensors
plot(fcdata_samples{1,2},diffdraw_fcleft_av, '-b')
grid on % add labels 
xlabel('Time [s]')
ylabel('Field intensity [in uV]')
title('Averaged centro-frontal left sensors for 0.6 probability condition')
hold on % to also add cond b
plot(fcdata_samples{1,2}, diffurn_fcleft_av,'-r','MarkerSize',6),
legend('0.6-draw  ','0.6-urn  ','Location','northwest') % Add a legend in the upper left:

%%% plot averaged parietal ERP for 0.6 right sensors
plot(fcdata_samples{1,2},diffdraw_fcright_av, '-b')
grid on % add labels 
xlabel('Time [s]')
ylabel('Field intensity [in uV]')
title('Averaged fronto-central right sensors for 0.6 probability condition')
hold on % to also add cond b
plot(fcdata_samples{1,2}, diffurn_fcright_av,'-r','MarkerSize',6),
legend('0.6-draw  ','0.6-urn  ','Location','northwest') % Add a legend in the upper left:

%%%% PLOT LEFT & RIGHT PARIETAL SENSORS for 0.8 condition
% First seperate trials/conditions
easydraw_fcleft     = data(1:4,:,3); % 1. easy-draw left trial
easydraw_fcright    = data(5:7,:,3); % 2. easy-draw right trial
easyurn_fcleft      = data(1:4,:,4); % 3. wasy-urn left trial
easyurn_fcright     = data(5:7,:,4); % 4. easy-urn right trial

% average over sensors
easydraw_fcleft_av  = mean(easydraw_fcleft,1);
easyurn_fcleft_av   = mean(easyurn_fcleft,1);
easydraw_fcright_av = mean(easydraw_fcright,1);
easyurn_fcright_av  = mean(easyurn_fcright,1);

%%% plot averaged parietal ERP for 0.8 left sensors
plot(fcdata_samples{1,2},easydraw_fcleft_av, '-b')
grid on % add labels 
xlabel('Time [s]')
ylabel('Field intensity [in uV]')
title('Averaged fronto-central left sensors for 0.8 probability condition')
hold on % to also add cond b
plot(fcdata_samples{1,2}, easyurn_fcleft_av,'-r','MarkerSize',6),
legend('0.8-draw  ','0.8-urn  ','Location','northwest') % Add a legend in the upper left:

%%% plot averaged parietal ERP for 0.8 right sensors
plot(fcdata_samples{1,2},easydraw_fcright_av, '-b')
grid on % add labels 
xlabel('Time [s]')
ylabel('Field intensity [in uV]')
title('Averaged fronto-central right sensors for 0.8 probability condition')
hold on % to also add cond b
plot(fcdata_samples{1,2}, easyurn_fcright_av,'-r','MarkerSize',6),
legend('0.8-draw  ','0.8-urn  ','Location','northwest') % Add a legend in the upper left:


