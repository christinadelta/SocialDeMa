%% PREPROCESSING & ANALYSIS OF ECONOMIC AND FACIAL ATTRACTIVENESS (BEST-CHOICE) EEG DATA

% VERSION 1 of formal analysis plan of the best-choice tasks EEG data with SPM12 
% Preprocessing/analysis script was created in February 2023.

% BEST-CHOICE TASKS:
% 1. economic (best-choice) task
% 2. facial attractiveness (best-choice) task

% Details of Preprocessing steps can be found in the doc file:
% https://docs.google.com/document/d/1JAkRF-O4qCPNAXVtlYUEgGkUumYOu4-yJQ0xP9ycKbM/edit?usp=sharing
% 
% The EEG data for the best-choice tasks were recorded in two .bdf files -- one file
% for each block.
% for preprocessing, up to epoching, each file is pre-processed seperately 

%%%%% VERY IMPORTANT %%%%%: 
% there were cases that a participant would have a noisy electrode. Quality
% of the EEG data during recording was logged in this file:
% https://docs.google.com/spreadsheets/d/16vhhOgA19vZDzn-K2K8PTeD3FnRPTRPAypzS450-NSk/edit#gid=0

% all preprocessing steps and analyses are done using this script. Meaning that we call all
% the SPM functions using this script.
% There shouldn't be a need to run any of the functions using the SPM eeg
% GUI. 

% The are three functions that are modified or created specifically for the
% best-choice tasks in the 'utilities' dir. These are:

% 1. createMontage.m -- This function creates a montage (performs channel
% selection) and re-refernces the signal to the average of all electrodes.
% In pre-processing we reference the data to the average of all electrodes.
% This is one of the most common and well known methods. This method
% requires excluding noisy electrodes from referencing. Use this log file
% to see if a given participant has a bad channel. In case there is a bad
% channel exclude it from the createMontage.m function (for more info see
% the documentation within the createMontage.m file).

% 2. bcTrialdef.m -- this function re-writes the events as as sample
% choices (reject choices) and accept choices. During the
% recording there was no way to code the triggers as reject and accept choices
% (because the trigger is sent right after the bead presentation screen is
% flipped), thus we need to do this change "manually" before epoching. For
% more info on how this is accomplished, see the documentation of the 
% bcTrialdef.m file. 

% 3. bcEstimateContrast.m -- this function is used to estimate contrasts
% for ERP and TF-power meeg objects.

%% Create required directories and define paths %%

% Add SPM12 to the matlab path if needed:
% addpath /Users/christinadelta/software/neuroscience/spm12 % change this to your own SPM12 path
% savepath

% create spm directories that will be used for storing the output MEEG
% objects mat files and images. MEEG objects are stored in the output dir,
% .mat files in the jobs dir and images are stored in directories that are
% created automatically during conversion of the MEEG objects to .nii
% files.

% clear workspace
clear all
clc

basedir     = pwd;

% if output and jobs directories do not exist, create them:
% outDir      = fullfile(basedir, 'output');
% jobsDir     = fullfile(basedir, 'jobs');
% cropDir     = fullfile(basedir, 'cropped');
outDir      = '/Volumes/DeepSpaceStuff/optimal_stopping_data/face/spm_analysis/output';
jobsDir     = '/Volumes/DeepSpaceStuff/optimal_stopping_data/face/spm_analysis/jobs';
cropDir     = '/Volumes/DeepSpaceStuff/optimal_stopping_data/face/spm_analysis/cropped';

if ~exist(outDir, 'dir') && ~exist(jobsDir, 'dir')
    mkdir(outDir)
    mkdir(jobsDir)
    mkdir(cropDir)
end

% add directories to the path
addpath(genpath(jobsDir));
addpath(genpath(fullfile(basedir, 'scripts')));
addpath(genpath(fullfile(basedir, 'utilities')));

% set data path
datadir         = '/Volumes/DeepSpaceStuff/optimal_stopping_data/data/facial_attractiveness/eeg/';  % change to your own data path 
subs            = dir(fullfile(datadir, '*sub*'));                      % how many subject folders?
nsubs           = length(subs);                                         % how many subjects?

taskname        = 'face';
blocks          = 2;
% conditions      = 2;
choices         = 2;
analysest       = 2; % analyses types (erp & tf)
nconstrasts     = 3;
condtypes       = {'reject' 'accept'};
contrastpref    = {'wra_' 'wr_' 'wa_'};


%% Preprocessing steps [ERPs: 1 - 12, TFRs: 1 - 15] %%

% these steps are used to preprocess the raw .bdf files, for each
% participant. Given that we will run evoked analysis and time-frequency
% analysis, we will low-pass filter the data in two ways, one for evoked
% analysis and one for TF analysis. Both MEEG objects will be preprocessed
% up until merging, however, only the ERPs MEEG object will undergo
% artefact rejection (the last step of preprocessing).

% init spm 
spm('defaults', 'eeg');

% loop over subjects
for sub = 1:nsubs


    fprintf('loading face block data\n')  
    subject         = subs(sub).name;
    subdir          = fullfile(datadir,subject);
    fprintf('\t reading data from subject %d\n',sub); 
    
    % create a subject sub-directory in outerps & outtfrs to store
    % subjected specific MEEG objects
    subout          = fullfile(outDir, sprintf('sub-%02d', sub));
    subjobs         = fullfile(jobsDir, sprintf('sub-%02d', sub));
    
    if ~exist(subout, 'dir') && ~exist(subjobs, 'dir')
        mkdir(subout)
        mkdir(subjobs)
    end

    % loop over blocks
    for block = 1:blocks

        fprintf('\t\t loading block %d\n\n',block);
        blockfile       = fullfile(subdir, sprintf('sub_%02d_%s_block_%02d.bdf', sub, taskname, block));

        %% STEP 1. convert the bdf file to MEEG object
        % create S struct for conversion
        S                   = [];
        S.dataset           = blockfile;
        S.mode              = 'continuous';
        S.channels          = 'all';

        % convert bdf file to spm object and D struct
        S.eventpadding      = 0;
        S.blocksize         = 3276800;
        S.checkboundary     = 1;
        S.saveorigheader    = 0;
        S.outpath           = fullfile(subout, sprintf('spmeeg_sub_%02d_%s_block_%02d.mat', sub, taskname, block));
        S.outfile           = S.outpath;
        S.timewin           = [];
        S.conditionlabels   = {'Undefined'};
        S.inputformat       = [];
        D                   = spm_eeg_convert(S); % convert raw data

        %% STEP 2. Create montage & re-reference 
       
        % create S struct for 
        S                   = [];
        S.D                 = fullfile(subout, sprintf('spmeeg_sub_%02d_%s_block_%02d.mat', sub, taskname, block));
        S.jobpath           = subjobs;
        S.block             = block;
        
        % Run createMontageFacialAttract.m function.
        % This function re-references by averaging across all electrodes. However,
        % an initial step requires knowing in advance if there is a noisy channel
        % and exclude it from averaging. Montage creation thus, need to be done
        % individually for every subject (e.g. pilot sub has one noisy channel [channel
        % 25]. This needs to be removed from averaging when re-referencing.
        % The createMontage.m file is in basedir/utilities; the function
        % needs to be modified manually 
        S                   = createMontageFacialAttract(S);
        
        S.mode              = 'write';
        S.blocksize         = 655360;
        S.prefix            = 'M';
        S.montage           = fullfile(S.jobpath, 'montage.mat');
        S.keepothers        = 0;
        S.keepsensors       = 1;
        S.updatehistory     = 1;
        D                   = spm_eeg_montage(S);

        %% STEP 3. High-pass filter 
        % Init S struct
        S                   = [];
        S.D                 = fullfile(subout, sprintf('Mspmeeg_sub_%02d_%s_block_%02d.mat', sub, taskname, block));
        S.type              = 'butterworth';
        S.band              = 'high';
        S.freq              = 0.1;
        S.dir               = 'twopass';
        S.order             = 5;
        S.prefix            = 'f';
        D                   = spm_eeg_filter(S);
        
        %% STEP 4. Downsample
        % Init S struct
        S                   = [];
        S.D                 = fullfile(subout, sprintf('fMspmeeg_sub_%02d_%s_block_%02d.mat', sub, taskname, block));
        S.fsample_new       = 256;
        S.method            = 'resample';
        S.prefix            = 'd';
        D                   = spm_eeg_downsample(S);


        %% STEP 5. Low-pass filter 
        
        % at this stage we will low-pass filter the downsampled file twice:
        % 1. for ERP analysis with frequency cutoff: 30Hz
        % 2. for TF analysis with frequency cutoff: 110Hz
        
        for i = 1:analysest
            S               = [];
            S.D             = fullfile(subout, sprintf('dfMspmeeg_sub_%02d_%s_block_%02d.mat', sub, taskname, block));
            S.type          = 'butterworth';
            S.band          = 'low';
            S.dir           = 'twopass';
            S.order         = 5;
            
            if i == 1 % if it's for ERP analysis
                S.freq      = 30;
                S.prefix    = 'erpf';
                
            else % if it's tf analysis
                S.freq      = 110;
                S.prefix    = 'tfrf';
            end
            
            D               = spm_eeg_filter(S);
            
        end % end of analysis types loop
        
        %% STEP 6. Epoch data 
       
        % FIRST SPECIFY TRIALS -- CREATE TRIAL DEFINITION MAT FILE
        % this part creates a mat file (trial definition) using the the previously
        % saved MEEG object (low-pass filtered). This mat file should contain:
        % 1. source of data
        % 2. time window
        % 3. trialdef (condition labels, event types, event values, trl shift)
        % 4. condition labels (cell with conditions/strings)
        % 5. trl (trial start, trial end, offset)
        
        % In the previous step we create two MEEG objects (one for ERP and
        % one for TF analysis). Epoching will be performed to each of the
        % two objects seperately as the tf-specific object does not require
        % baseline-correction. 

        for i = 1:analysest

            % init S struct
            S                               = [];
            S.trialdef(1).conditionlabel    = 'reject';
            S.trialdef(1).eventtype         = 'STATUS';
            S.trialdef(1).eventvalue        = 1;
            S.trialdef(1).trlshift          = 0;
            S.trialdef(2).conditionlabel    = 'accept';
            S.trialdef(2).eventtype         = 'STATUS';
            S.trialdef(2).eventvalue        = 2;
            S.trialdef(2).trlshift          = 0;
            S.timewin                       = [-500 800];
            S.eventpadding                  = 0;
            S. prefix                       = 'e';
            if i == 1 % if this is the ERP object
                
                S.D                         = fullfile(subout, sprintf('erpfdfMspmeeg_sub_%02d_%s_block_%02d.mat', sub, taskname, block));
                S.bc                        = 1;
                
            else % if this is the TFR object
                
                S.D                         = fullfile(subout, sprintf('tfrfdfMspmeeg_sub_%02d_%s_block_%02d.mat', sub, taskname, block));
                S.bc                        = 0;
            end

            % run the trial definition function here to get the trl matrix and
            % condition labels
            [trl, conditionlabels, S]       = faTrialdef(S, sub);
            
            D                               = spm_eeg_epochs(S);

        end % end of analysis set

        % This is the last preprocessing step at the block level. Now the
        % files need to be merged and block files can be deleted since they
        % are not going to be used again. 
        

    end % end of blocks loop

    %% STEP 7. Merge all block objects 
    
    for i = 1:analysest
        
        % init S struct
        S                   = [];
        S.recode.file       = '.*';
        S.recode.labelorg   = '.*';
        S.recode.labelnew   = '#labelorg#';
        S.prefix            = 'c';
%         S.save              = 1;
        
        if i == 1 % if this is the ERP object
            S.D             = [fullfile(subout, sprintf('eerpfdfMspmeeg_sub_%02d_face_block_01.mat', sub))
                fullfile(subout, sprintf('eerpfdfMspmeeg_sub_%02d_face_block_02.mat', sub))
                ];
        else
            S.D             = [fullfile(subout, sprintf('etfrfdfMspmeeg_sub_%02d_face_block_01.mat', sub))
                fullfile(subout, sprintf('etfrfdfMspmeeg_sub_%02d_face_block_02.mat', sub))
                ];
        end
        
        D                   = spm_eeg_merge(S);
        
    end % end of analysis type loop

    % move .mat & .dat files that were created and stored in pwd to the
    % subout directory
    movefile('*.mat', subout); movefile('*.dat', subout);

    %% STEP 8. Define cordinates/channel locations 
%     
%     % this step is not 100% necessary as default biosemi cordinates are
%     % used during montage creation. I include it just to be sure..
%     % but it can be commented-out
% %     for i = 1:analysest
% %         
% %         % prepare defauld spm cordinates/locations
% %         S               = [];
% %         S.task          = 'defaulteegsens';
% %         S.save          = 1;
% %         
% %         if i == 1 % if this is the ERP object
% %             
% %             S.D         = fullfile(subout, sprintf('ceerpfdfMspmeeg_sub_%02d_economic_block_01.mat', sub));
% %  
% %         else
% %             S.D         = fullfile(subout, sprintf('cetfrfdfMspmeeg_sub_%02d_economic_block_01.mat', sub));
% % 
% %         end
% %         
% %         D               = spm_eeg_prep(S); % run eeg prep function
% %         
% %         % set coordinates 
% %         S               = [];
% %         S.task          = 'setcoor2d';
% %         
% %         if i == 1 % if this is the ERP object
% %             
% %             S.D         = fullfile(subout, sprintf('ceerpfdfMspmeeg_sub_%02d_economic_block_01.mat', sub));
% %  
% %         else
% %             S.D         = fullfile(subout, sprintf('cetfrfdfMspmeeg_sub_%02d_economic_block_01.mat', sub));
% % 
% %         end
% %         
% %         S.xy            = [0.358127431070432 0.230361291206221 0.364687127498288 0.412094935710624 0.319499558310764 0.220586712887751 0.121753409066311 0.0551811363296929 0.178452330461275 0.296722753529597 0.402418214252759 0.402232139889891 0.295640676509763 0.176774133670791 0.05 0.0952050213281216 0.204129030986282 0.31059947251063 0.407784984615512 0.419107072168501 0.337261391153001 0.255936320417754 0.176983444035271 0.0990531376307864 0.274933610032685 0.374805341575767 0.381476047145341 0.485128199634413 0.490099574981175 0.493654015274171 0.496027125746896 0.497486874012439 0.494127551616874 0.63115333565553 0.764761267196509 0.633966624967236 0.497096812715531 0.498387742743428 0.590119958154716 0.684446454131226 0.784927979314371 0.878951869887076 0.947194362550941 0.826079267890216 0.706962412478808 0.597352645199803 0.498716216854912 0.498366895939448 0.598675341202977 0.70801035314695 0.826759976212305 0.95 0.898334363412286 0.795911344188124 0.69078343536388 0.593172800977718 0.57877084537302 0.658765804468392 0.735389719899426 0.810437806844084 0.879007062596769 0.708438648401535 0.611863922610018 0.599476237529535
% %         0.94291962707719 0.889006557470294 0.844625065001036 0.711929958683174 0.727189162575643 0.752686038230401 0.785972668551047 0.638824390760336 0.622237843107835 0.612055755800727 0.604719623778926 0.504717514961853 0.499865625106497 0.489925660632489 0.477921883184179 0.346661648796377 0.375837233637139 0.399737181066574 0.41403150193654 0.321035174889695 0.30567474167612 0.277954382109169 0.241048698933626 0.191392522403548 0.173810943732217 0.218490107582391 0.139058031510229 0.05 0.134825637229592 0.229144765608567 0.325406282673219 0.415924315197709 0.95 0.947839059007246 0.895554940224051 0.844643723365574 0.824931733550042 0.708258571334129 0.714599145507353 0.73251653002839 0.754709519350322 0.7891126332849 0.637244898999659 0.620328254919284 0.610930378565758 0.604973846888264 0.602533725218943 0.506171957423089 0.501865455408111 0.495100231864275 0.484556343615868 0.468273841263323 0.33147104080173 0.369228471090915 0.394026914867007 0.410294741704492 0.320570318800589 0.301276182177249 0.266698929482614 0.226857688879681 0.170533558640877 0.16331740311369 0.210477257752421 0.133428670630485];
% %         
% %         S.label = {'Fp1', 'AF7', 'AF3', 'F1', 'F3', 'F5', 'F7', 'FT7', 'FC5', 'FC3', 'FC1', 'C1', 'C3', 'C5', 'T7',...
% %             'TP7', 'CP5', 'CP3', 'CP1', 'P1', 'P3', 'P5', 'P7', 'P9', 'PO7', 'PO3', 'O1', 'Iz', 'Oz', 'POz', 'Pz', 'CPz',...
% %             'Fpz', 'Fp2', 'AF8', 'AF4', 'AFz', 'Fz', 'F2', 'F4', 'F6', 'F8', 'FT8', 'FC6', 'FC4', 'FC2', 'FCz','Cz', 'C2',...
% %             'C4', 'C6', 'T8', 'TP8', 'CP6', 'CP4', 'CP2', 'P2', 'P4', 'P6', 'P8', 'P10', 'PO8',  'PO4', 'O2'};
% %         
% %         %save the file 
% %         S.save          = 1;
% %         D               = spm_eeg_prep(S);
% %         
% %     end
%  

    %% STEP 9: Artefact rejection
    
    % artefact rejection is performed only on the ERP-specific MEEG object.
    % 
    S                               = [];
    S.D                             = fullfile(subout, sprintf('ceerpfdfMspmeeg_sub_%02d_face_block_01.mat', sub));
    S.mode                          = 'reject';
    S.badchanthresh                 = 0.2;
    S.prefix                        = 'a';
    S.append                        = true;
    S.methods.channels              = {'EEG'};
    S.methods.fun                   = 'threshchan';
    S.methods.settings.threshold    = 120;
    S.methods.settings.excwin       = 1000;
    D                               = spm_eeg_artefact(S);

    %% STEP 10: Average conditions (Averaged ERPs)
    
    % condition averageing is performed only on the ERP-specific MEEG object.
    S                       = [];
    S.D                     = fullfile(subout, sprintf('aceerpfdfMspmeeg_sub_%02d_face_block_01.mat', sub));
    S.robust.ks             = 3;
    S.robust.bycondition    = true;
    S.robust.savew          = false;
    S.robust.removebad      = false;
    S.circularise           = false;
    S.prefix                = 'm';
    D                       = spm_eeg_average(S);

    %% Step 10 Optional: Re-apply low pass filter (optional)
    
    S               = [];
    S.D             = fullfile(subout, sprintf('maceerpfdfMspmeeg_sub_%02d_face_block_01.mat', sub));
    S.type          = 'butterworth';
    S.band          = 'low';
    S.dir           = 'twopass';
    S.order         = 5;
    S.freq          = 30;
    S.prefix        = '';
    D               = spm_eeg_filter(S);

    %% STEP 11: contrast ERP averaged conditions
    
    % We will first contrast averaged ERPs. We will create 3 contrasts as:
    % accpet vs reject, only reject, only accept 
    S                       = [];
    S.D                     = fullfile(subout, sprintf('maceerpfdfMspmeeg_sub_%02d_face_block_01.mat', sub));
    S.object                = spm_eeg_load(S.D);
    tmp                     = S.object.conditions;
    
    % estimate ERP contrasts
    faEstimateContrasts(S, tmp, sub)


    %% STEP 12: Convert to 3D volumes 

    % convert the 3 contrast MEEG objects into .nii file
    for c = 1:nconstrasts
        S               = [];
        S.D             = fullfile(subout, sprintf('%smaceerpfdfMspmeeg_sub_%02d_face_block_01.mat', contrastpref{c}, sub));
        S.mode          = 'scalp x time';
        S.conditions    = {};
        S.channels      = 'EEG';
        S.timewin       = [-Inf Inf];
        S.freqwin       = [-Inf Inf];
        S.prefix        = '';
    
        D               = spm_eeg_convert2images(S);
    end % end of contrasts loop

    %% RUN TIME-FREQUENCY ANALYSIS %%

%     at this step we will start using the TF-specific MEEG file that
%     hasn't been used after the merging step. A lot of the
%     preprocessing steps for TF representations will be similar to ERPs
%     (e.g., condition averaging, contrast creation and convertion to
%     3D images), however, the first 3 TF-specific steps will be
%     different.
%     Here we continue with step 9, since definition of channel locations was step 8.

%% STEP 9: Time-Frequency Morlet Decomposition

%     at this step SPM12 creates 2 objects, one for power (tf) and one
%     for phase (ph). At the remaining of the TF preprocessing steps we will only use
%     the tf file. The ph file can be deleted.

    S                       = [];
    S.D                     = fullfile(subout, sprintf('cetfrfdfMspmeeg_sub_%02d_face_block_01.mat', sub));
    S.channels              = {'all'};
    S.frequencies           = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55];
    S.timewin               = [-500 800];
    S.phase                 = 1;
    S.method                = 'morlet';
    S.settings.ncycles      = 7;
    S.settings.timeres      = 0;
    S.settings.subsample    = 2; % we downsampled by a factor of 2
    S.prefix                = '';
    D                       = spm_eeg_tf(S);

    %% STEP 10: Baseline rescaling (only power file) 

    S                       = [];
    S.D                     = fullfile(subout, sprintf('tf_cetfrfdfMspmeeg_sub_%02d_face_block_01.mat', sub));
    S.method                = 'LogR';
    S.prefix                = 'r';
    S.timewin               = [-500 -50];
    S.pooledbaseline        = 0;
    D                       = spm_eeg_tf_rescale(S);

    %% STEP 11: Average power over frequency

    % before averaging conditions, we will need to average the TF
    % object over beta frequency (as this is the focus of the analysis)

    S                       = [];
    S.D                     = fullfile(subout, sprintf('rtf_cetfrfdfMspmeeg_sub_%02d_face_block_01.mat', sub));
    S.freqwin               = 13:30;
    S.prefix                = 'P';
    D                       = spm_eeg_avgfreq(S);

%       %% STEP 12: Average power over time
% 
%     % before averaging conditions, we will need to average the TF
%     % object over the entire peristimulus time (-500 to 800)
%     S                       = [];
%     S.D                     = fullfile(subout, sprintf('rtf_cetfrfdfMspmeeg_sub_%02d_face_block_01.mat', sub));
%     S.timewin               = [-500 800];
%     S.prefix                = 'S';
%     D                       = spm_eeg_avgtime(S);
    %% STEP 12: Average conditions

    S                       = [];
    S.D                     = fullfile(subout, sprintf('Prtf_cetfrfdfMspmeeg_sub_%02d_face_block_01.mat', sub));
    S.robust.ks             = 3;
    S.robust.bycondition    = true;
    S.robust.savew          = false;
    S.robust.removebad      = false;
    S.circularise           = false;
    S.prefix                = 'm';
    D                       = spm_eeg_average(S);

    %% STEP 13: Compute contrasts of averaged power file

    % As with evoked responses, TFRs also will be contrasted as: accept vs
    % reject, only accept, only reject
    S                       = [];
    S.D                     = fullfile(subout, sprintf('mPrtf_cetfrfdfMspmeeg_sub_%02d_face_block_01.mat', sub));
    S.object                = spm_eeg_load(S.D);
    tmp                     = S.object.conditions;
    
    % estimate ERP contrasts
    faEstimateContrasts(S, tmp, sub)
    
    %% STEP 14: Convert to 3D volumes 

    % convert the 3 contrast MEEG objects into .nii file
    for c = 1:nconstrasts
        S               = [];
        S.D             = fullfile(subout, sprintf('%smPrtf_cetfrfdfMspmeeg_sub_%02d_face_block_01.mat', contrastpref{c}, sub));
        S.mode          = 'scalp x time';
        S.conditions    = {};
        S.channels      = 'EEG';
        S.timewin       = [-Inf Inf];
        S.freqwin       = [-Inf Inf];
        S.prefix        = '';
    
        D               = spm_eeg_convert2images(S);

    end

end % end of subjects loop

% THIS IS IT WITH PREPROCESSING!!!

%% Compute grand averages

