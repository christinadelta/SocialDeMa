%% BEADS EEG DATA - FORMAL PREPROCESSING SCRIPT W/ FIELDTRIP

% DATA are only preprocessed in FieldTrip, we save the resukting ERPs and
% TFRs in structures for statistical analyeses in SPM 

% christinadelta 
% Created: Jan 2022
% Version2: 13/03/2022
% Version3: 18/05/2022

% final preprocessing of the EEG data - beads task using FIELDTRIP 

% A few differences in pre-processing methods between biosemi (.bdf) data
% and other devises is that:
% 1) with biosemi we don't need to do re-referencing to obtain EOG channels,
% the signal is recorded as bipolar 
% 2) biosemi records and saves the signal UN-REFERENCED, thus, we don't
% include an "implicitref" during preprocessing. We only need to provide the "reref" channel(s)

% HOW DO WE COMPUTE ERPs (and TFRs):
% In fieldtrip epochs are called "trials"! In our information sampling
% (beads) task, every sequence can have up to 10
% samples/epochs/trials(fieltrip).

% OUR MAIN CONDITIONS ARE: 
% A) COLOUR PROBABILITY: EASY/0.8 & DIFFICULT/0.6
% B) CHOICE TYPES: DRAW CHOICES & URN CHOICES

% PREPROCESSING PIPELINE:
% load continues EEG data - keep only relevant trials in seperate
% preprocessing structures
% Re-reference the data
% Run ICA and remove components containing EOG artefacts
% Visualise all channels and interpolate (if needed)
% Downsample to 256Hz
% merge the different block/run trials into one struct 
% filter & epoch trials (this is called timelock analysis in FT). We will extract 3 types of epochs:
% 1) separate epochs for difficult and easy conditions (1 epoch per
% stimulus trigger)
% 2) seperate epochs per choice types (1 epoch per choice type: draw - urn choice)
% 3) combined conditions epochs: easy/draw, easy/urn, diff/draw, diff/urn

% WILL PROBABLY SAVE EPOCHED DATA AND RUN ERPs & TFRs in SPM

%% Define paths and variables 

% change paths as appropriate:
basedir             = '/Users/christinadelta/Desktop/os_data';
task                = 'beads';
taskdir             = fullfile(basedir, task); 

subs                = dir(fullfile(taskdir, '*sub*'));
nsubs               = length(subs);

subname             = {subs.name}; % only keep the subject name 

% define initial variables 
blocks              = 4; 
blocktrials         = 13; 
totaltrials         = blocks * blocktrials;
totalconds          = 4; % easy/draw, easy/urn, diff/draw, diff/urn
nconds              = 2; % [easy, diff] and [draw choices, urn choices]
conditions          = {'Easy', 'Difficult'};

% define layout-montage for topoplots (this will be mainly used for ICA
% plots
% check the biosemi64 layout that will be used for plotting topomaps
cfglayout           = [];
cfglayout.layout    = 'biosemi64.lay';
lo                  = ft_layoutplot(cfglayout);

%%  Load the data 

% loop over subjects 
for subI = 1:nsubs 
    
    % specify subI dir 
    subIdir = fullfile(taskdir, subname{subI});
    
    % loop over blocks 
    for blockI = 1:blocks
        
        % load subI blockI .bdf file 
        subFile = fullfile(subIdir, sprintf('sub_%02d_%s_block_%02d.bdf', subI, task, blockI));
        cfg                     = [];
        cfg.dataset             = subFile;

        cfg.trialdef.eventtype  = 'STATUS';
        cfg.trialdef.eventvalue = [1 2 3 4 102 103];
        cfg.trialdef.prestim    = 0.5;
        cfg.trialdef.poststim   = 0.8; 
        cfg                     = ft_definetrial(cfg);
        
        %% Split block-trials in different preprocessing data structures (useful for timelock analysis/epoching)
        
        % extract trial list
        trl                     = cfg.trl;
        tstart                  = 102; % start of sequence
        tend                    = 103; % end of sequence

        blocksequences          = length(find(trl(:,4) == tstart));
        trialstart              = find(trl(:,4) == tstart); % should be 13 
        trialend                = find(trl(:,4) == tend);   % should be 13 

        counter                 = 0;
        
        % loop over block sequences 
        for iTrial = 1:blocksequences
            
            % get sequence(iTrial) length (this stores the number of bead presentations in iTrial sequence) 
            tmp                 = length(trialstart(iTrial)+1: trialend(iTrial)-1); 
            
            % store trial/sequence and block numbers
            for j = 1:tmp
                
                cnt             = counter + j; 
                trialnum(cnt,1) = iTrial;
                trialnum(cnt,2) = blockI;
            end 
            
            % update counter 
            counter             = counter + tmp ;
            
            % clear workspace from not used vars
            clear tmp cnt
            
            
        end % end of trials loop
        
        % remove trialstart and trialend from the trl list
        trl(trl(:,4) == tstart, :)  = [];
        trl(trl(:,4) == tend, :)    = [];
        trl(:,5)                    = trialnum(:,1); % add trialnum to the main list
        trl(:,6)                    = trialnum(:,2); % add trialnum to the main list
        
        % now split data into conditions for later preprocessing, & move condition data to
        % new data structures     
        trl_length                      = length(trl);
        
        % add a condition column (1: easy, 2: diff)
        for i = 1:trl_length

            if trl(i,4) == 1 | trl(i,4) == 2
                trl(i,7) = 1;

            elseif trl(i,4) == 3 | trl(i,4) == 4
                trl(i,7) = 2;

            end
        end
        
        % clear vars for memory eficiency
        clear trialend trialstart j i counter cnt trialnum 
        
        % split trial data in easy and difficult conditions 
        trl_easy            = trl((find(trl(:,7) == 1)),:);
        trl_diff            = trl((find(trl(:,7) == 2)),:);
        
        % split trl data into choice condtions [draw choice - urn choice] 
        draw                = 0;
        urn                 = 0;
        c                   = 0; % counter index
        l                   = 1; % last draw index (urn choice)
        
        for icond = 1:nconds
            for itrial = 1:blocksequences

                tmp = find(trl(:,7)== icond & trl(:,5)== itrial);

                if ~isempty(tmp)
                    tl                  = length(tmp)-1;
                    draw(c+1:c+tl)      = tmp(1:end-1); % only pick 
                    urn(:,l)            = tmp(end);

                    % update c and l 
                    c                   = c + tl;
                    l                   = l + 1;
                    
                end % end of if statement
            end % end of trials loop
        end % end of iconds loop
        
        trl_draw                        = trl((draw),:);
        trl_urn                         = trl((urn),:);
        
        % further split draw and urn trls into easy and diff draw and
        % urn
        tmp_draweasy                    = find(trl_draw(:,7) == 1);
        tmp_drawdiff                    = find(trl_draw(:,7) == 2);
        trl_draweasy                    = trl_draw((tmp_draweasy),:);
        trl_drawdiff                    = trl_draw((tmp_drawdiff),:);
        
        tmp_urneasy                     = find(trl_urn(:,7) == 1);
        tmp_urndiff                     = find(trl_urn(:,7) == 2);
        trl_urneasy                     = trl_urn((tmp_urneasy),:);
        trl_urndiff                     = trl_urn((tmp_urndiff),:);
        
        clear draw urn tmp itrial icondc l tl tmp_urneasy tmp_urndiff tmp_draweasy tmp_drawdiff
        
        %% re-reference/preprocess/resample

        cfg.reref           = 'yes';
        cfg.refchannel      = {'EXG1' 'EXG2'};
        cfg.demean          = 'yes';
%         cfg.hpfilter        = 'yes';        % enable high-pass filtering
%         cfg.lpfilter        = 'yes';        % enable low-pass filtering
%         cfg.hpfreq          = 20;           % set up the frequency for high-pass filter
%         cfg.lpfreq          = 250;          % set up the frequency for low-pass filter
%         cfg.bpfilter        = 'yes'; % enable bandpass filtering
%         cfg.bpfreq          = [0.16 100];
        cfg.baselinewindow  = [-0.500 -0.050];
        
        % re-write the trl list to the cfg struct and preprocess all data, easy/diff data & draw/urn choice data structs 
        % the only thing that is changing here is "cfg.trl"
        cfg.trl             = trl;
        alldata             = ft_preprocessing(cfg);
        
        % resample alldata to 256Hz
        cfg_res             = [];
        cfg_res.resamplefs  = 256;
        alldata             = ft_resampledata(cfg_res, alldata);
        
        % preprocess easy data
        cfg.trl             = trl_easy;
        easydata            = ft_preprocessing(cfg);
        
        % resample easy data to 256Hz
        cfg_res             = [];
        cfg_res.resamplefs  = 256;
        easydata            = ft_resampledata(cfg_res, easydata);
        
        % preprocess diff data
        cfg.trl             = trl_diff;
        diffdata            = ft_preprocessing(cfg);
        
        % resample diff data to 256Hz
        cfg_res             = [];
        cfg_res.resamplefs  = 256;
        diffdata            = ft_resampledata(cfg_res, diffdata);
        
        % preprocess draw choice data 
        cfg.trl             = trl_draw;
        drawdata            = ft_preprocessing(cfg);
        
        % resample to 256Hz
        cfg_res             = [];
        cfg_res.resamplefs  = 256;
        drawdata            = ft_resampledata(cfg_res, drawdata);
        
        % preprocess urn choice data
        cfg.trl             = trl_urn;
        urndata             = ft_preprocessing(cfg);
        
        % resample to 256Hz
        cfg_res             = [];
        cfg_res.resamplefs  = 256;
        urndata             = ft_resampledata(cfg_res, urndata);
        
        % preprocess easy draw choice data 
        cfg.trl             = trl_draweasy;
        draweasydata        = ft_preprocessing(cfg);
        
        % resample to 256Hz
        cfg_res             = [];
        cfg_res.resamplefs  = 256;
        draweasydata        = ft_resampledata(cfg_res, draweasydata);
        
        % preprocess diff draw choice data 
        cfg.trl             = trl_drawdiff;
        drawdiffdata        = ft_preprocessing(cfg);
        
        % resample to 256Hz
        cfg_res             = [];
        cfg_res.resamplefs  = 256;
        drawdiffdata        = ft_resampledata(cfg_res, drawdiffdata);
        
        % preprocess easy urn choice data 
        cfg.trl             = trl_urneasy;
        urneasydata         = ft_preprocessing(cfg);
        
        % resample to 256Hz
        cfg_res             = [];
        cfg_res.resamplefs  = 256;
        urneasydata         = ft_resampledata(cfg_res, urneasydata);
        
        % preprocess diff urn choice data 
        cfg.trl             = trl_urndiff;
        urndiffdata         = ft_preprocessing(cfg);
        
        % resample to 256Hz
        cfg_res             = [];
        cfg_res.resamplefs  = 256;
        urndiffdata         = ft_resampledata(cfg_res, urndiffdata);
      
        
        % only keep eeg channels from now on % for all data structures
        % [exclude eog and mastoids]
        cfg                 = [];
        cfg.channel         = [1:64];
        
        alldata             = ft_selectdata(cfg, alldata);
        easydata            = ft_selectdata(cfg, easydata);
        diffdata            = ft_selectdata(cfg, diffdata);
        drawdata            = ft_selectdata(cfg, drawdata);
        urndata             = ft_selectdata(cfg, urndata);
        draweasydata        = ft_selectdata(cfg, draweasydata);
        drawdiffdata        = ft_selectdata(cfg, drawdiffdata);
        urneasydata         = ft_selectdata(cfg, urneasydata);
        urndiffdata         = ft_selectdata(cfg, urndiffdata);
        
        %% Store and save preprocessed block data
        
        % first store all data in data struct
        data.alldata        = alldata;
        data.easydata       = easydata;
        data.diffdata       = diffdata;
        data.drawdata       = drawdata;
        data.urndata        = urndata;
        data.draweasydata   = draweasydata;
        data.drawdiffdata   = drawdiffdata;
        data.urneasydata    = urneasydata;
        data.urndiffdata    = urndiffdata;
        
        % save preprocesssed block data in a .mat file 
        save(['beads_analysis/prepro/beads_preproc_sub_', num2str(subI), '_block_', num2str(blockI)], 'data')
        
        % clear workspace for memory efficiency 
        clear data alldata easydata diffdata drawdata urndata cfg trl_diff trl_easy trl_draw trl_urn
  
    end % end of block loop
    
    %% Append data into one structure
    
    % Load the matfiles and concatenate all blocks 
    for block = 1:blocks
        
        % if data structure is already in workspace, comment the part disp
        % and load parts 
        disp(['loading beads_analysis/prepro/beads_preproc_sub_', num2str(subI), '_block_', num2str(block)])
        load(['beads_analysis/prepro/beads_preproc_sub_', num2str(subI), '_block_', num2str(block)], 'data')
        
        partdata(block)         = data.alldata;
        easy_data(block)        = data.easydata;
        diff_data(block)        = data.diffdata;
        draw_data(block)        = data.drawdata;
        urn_data(block)         = data.urndata;
        draw_easydata(block)    = data.draweasydata;
        draw_diffdata(block)    = data.drawdiffdata;
        urn_easydata(block)     = data.urneasydata;
        urn_diffdata(block)     = data.urndiffdata;
        
          
    end
    
    % append blocks (merge)
    cfg             = [];
    alldata         = ft_appenddata(cfg, partdata(1), partdata(2), partdata(3), partdata(4));
    easydata        = ft_appenddata(cfg, easy_data(1), easy_data(2), easy_data(3), easy_data(4));
    diffdata        = ft_appenddata(cfg, diff_data(1), diff_data(2), diff_data(3), diff_data(4));
    drawdata        = ft_appenddata(cfg, draw_data(1), draw_data(2), draw_data(3), draw_data(4));
    urndata         = ft_appenddata(cfg, urn_data(1), urn_data(2), urn_data(3), urn_data(4));
    draweasydata    = ft_appenddata(cfg, draw_easydata(1), draw_easydata(2), draw_easydata(3), draw_easydata(4));
    drawdiffdata    = ft_appenddata(cfg, draw_diffdata(1), draw_diffdata(2), draw_diffdata(3), draw_diffdata(4));
    urneasydata     = ft_appenddata(cfg, urn_easydata(1), urn_easydata(2), urn_easydata(3), urn_easydata(4));
    urndiffdata     = ft_appenddata(cfg, urn_diffdata(1), urn_diffdata(2), urn_diffdata(3), urn_diffdata(4));
    
    % clear workspace
    clear partdata easy_data diff_data draw_data urn_data draw_easydata draw_diffdata urn_easydata urn_diffdata
    
    %% Run ICA algorithm for artefact removal
        
    % run ICA 
    cfg_ica                 = [];
    cfg_ica.method          = 'runica'; % which ica we'll run?
    cfg_ica.numcomponent    = 20;
    allcomp                 = ft_componentanalysis(cfg_ica, alldata);
    easycomp                = ft_componentanalysis(cfg_ica, easydata);
    diffcomp                = ft_componentanalysis(cfg_ica, diffdata);
    drawcomp                = ft_componentanalysis(cfg_ica, drawdata);
    urncomp                 = ft_componentanalysis(cfg_ica, urndata);
    draweasycomp            = ft_componentanalysis(cfg_ica, draweasydata);
    drawdiffcomp            = ft_componentanalysis(cfg_ica, drawdiffdata);
    urneasycomp             = ft_componentanalysis(cfg_ica, urneasydata);
    urndiffcomp             = ft_componentanalysis(cfg_ica, urndiffdata);


    %% Visualise components and remove if needed

    % specify number of components 
    ncomps = 1:20; % 20 componenets

    cfg=[];
    cfg.component   = ncomps; 
    cfg.layout      = 'biosemi64.lay';   
    cfg.zlim        = 'maxabs';

    % run the lines below one-by-one to get components for each dataset
    allica          = ft_topoplotIC(cfg, allcomp);
    easyica         = ft_topoplotIC(cfg, easycomp);
    diffica         = ft_topoplotIC(cfg, diffcomp);
    drawica         = ft_topoplotIC(cfg, drawcomp);
    urnica          = ft_topoplotIC(cfg, urncomp);
    draweasyica     = ft_topoplotIC(cfg, draweasycomp);
    drawdiffica     = ft_topoplotIC(cfg, drawdiffcomp);
    urneasyica      = ft_topoplotIC(cfg, urneasycomp);
    urndiffica      = ft_topoplotIC(cfg, urndiffcomp);

    % also visualise in time course (also visualise one-by-one)
%     cfg.viewmode    = 'vertical';
%     timecourse      = ft_databrowser(cfg, allcomp);
%     timecourse      = ft_databrowser(cfg, easycomp);
%     timecourse      = ft_databrowser(cfg, diffcomp);
%     timecourse      = ft_databrowser(cfg, drawcomp);
%     timecourse      = ft_databrowser(cfg, urncomp);

    % remove components that reflect eog artifacts
    cfg             = [];
    cfg.component   = [7 11 13]; % the exact numbers varies per run
    ica_alldata     = ft_rejectcomponent(cfg, allcomp);

    cfg             = [];
    cfg.component   = [9]; % the exact numbers varies per run
    ica_easydata    = ft_rejectcomponent(cfg, easycomp);

    cfg             = [];
    cfg.component   = [6]; % the exact numbers varies per run
    ica_diffdata    = ft_rejectcomponent(cfg, diffcomp);

    cfg             = [];
    cfg.component   = [8]; % the exact numbers varies per run
    ica_drawdata    = ft_rejectcomponent(cfg, drawcomp);

    cfg             = [];
    cfg.component   = [7 11 12]; % the exact numbers varies per run
    ica_urndata     = ft_rejectcomponent(cfg, urncomp);
    
    cfg                 = [];
    cfg.component       = [9]; % the exact numbers varies per run
    ica_draweasydata    = ft_rejectcomponent(cfg, draweasycomp);
    
    cfg                 = [];
    cfg.component       = [5]; % the exact numbers varies per run
    ica_drawdiffdata    = ft_rejectcomponent(cfg, drawdiffcomp);
    
    cfg             = [];
    cfg.component   = [8 11 15]; % the exact numbers varies per run
    ica_urneasydata = ft_rejectcomponent(cfg, urneasycomp);
    
    cfg             = [];
    cfg.component   = [7 10]; % the exact numbers varies per run
    ica_urndiffdata = ft_rejectcomponent(cfg, urndiffcomp);

    
    %% Visualise trials/channels and reject if needed 
        
    % if a channel needs to be rejected, should I also interpolate here?
    % visually inspect "channels" and exclude if needed -- all data 
    cfg                 = [];
    cfg.method          = 'channel';
    clean_alldata       = ft_rejectvisual(cfg, ica_alldata); % if ica components removed, use "ica_alldata"
    clean_easydata      = ft_rejectvisual(cfg, ica_easydata);
    clean_diffdata      = ft_rejectvisual(cfg, ica_diffdata);
    clean_drawdata      = ft_rejectvisual(cfg, ica_drawdata);
    clean_urndata       = ft_rejectvisual(cfg, ica_urndata);
    clean_draweasydata  = ft_rejectvisual(cfg, ica_draweasydata);
    clean_drawdiffdata  = ft_rejectvisual(cfg, ica_drawdiffdata);
    clean_urneasydata   = ft_rejectvisual(cfg, ica_urneasydata);
    clean_urndiffdata   = ft_rejectvisual(cfg, ica_urndiffdata);
    

    % visualise 
    cfg                 = [];
    ft_databrowser(cfg, clean_alldata)
    ft_databrowser(cfg, clean_easydata)
    ft_databrowser(cfg, clean_diffdata)
    ft_databrowser(cfg, clean_drawdata)
    ft_databrowser(cfg, clean_urndata)
    ft_databrowser(cfg, clean_draweasydata)
    ft_databrowser(cfg, clean_drawdiffdata)
    ft_databrowser(cfg, clean_urneasydata)
    ft_databrowser(cfg, clean_urndiffdata)

    %% Interpolate bad channels - maybe should do that after ft_preprocessing - will work on that 

%     % first specify neighbours
%     cfg_neighb              = [];
%     cfg_neighb.layout       = 'biosemi64.lay'; %in meters
%     cfg_neighb.method       = 'template'; % we could also use the template method 
%     neighbours              = ft_prepare_neighbours(cfg_neighb);
% 
%     % interpolate
%     cfg = [];
%     cfg.badchannel     = artif.badchannel;
%     cfg.method         = 'weighted';
%     cfg.neighbours     = neighbours;
%     inter_alldata      = ft_channelrepair(cfg,clean_alldata);
% 
%     % visualise interpolated data
%     cfg                 = [];
%     ft_databrowser(cfg, inter_alldata)

    %% Lowpass filter data and run timelock (ERPs) analysis (3 analyses) both averaged and with trials
    
    % AVERAGE TIMELOCK ANALYSES
    
    %% Save subject preprocessed/cleaned data
    
    % if sub folder doesn't exist, create one
    resultsfolder  = fullfile(basedir, 'beads_analysis', sprintf('sub-%02d', subI));
    
    if ~exist(resultsfolder, 'dir')
        mkdir(resultsfolder)
    end
    
    % save all epoched datasets for analyses in SPM
    cleandata.alldata       = clean_alldata;
    cleandata.easydata      = clean_easydata;
    cleandata.diffdata      = clean_diffdata;
    cleandata.drawdata      = clean_drawdata;
    cleandata.urndata       = clean_urndata;
    cleandata.draweasydata  = clean_draweasydata;
    cleandata.drawdiffdata  = clean_drawdiffdata;
    cleandata.urneasydata   = clean_urneasydata;
    cleandata.urndiffdata   = clean_urndiffdata;
    
    datasets                = fullfile(resultsfolder, 'cleandata');
    save(datasets, 'cleandata')
    
 
end % end of subject loop

