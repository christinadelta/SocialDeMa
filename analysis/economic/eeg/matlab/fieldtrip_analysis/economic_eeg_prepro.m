%% BEADS FORMAL PREPROCESSING AND ANALYSES SCRIPT W/ FIELDTRIP

% christinadelta 
% Date: Feb 2022
% last Update 27/02/2022

% final preprocessing and analyses of the EEG data - Best-Choice Economic task using
% FIELDTRIP 

% A few differences in pre-processing methods between biosemi (.bdf) data
% and other devises is that:
% 1) with biosemi with don't need to do re-referencing to obtain EOG channels,
% the signal is recorded as bipolar 
% 2) biosemi records and saves the signal UN-REFERENCED, thus, we don't
% include an "implicitref" during analysis. We only need to provide the "reref" channel(s)

% HOW DO WE COMPUTE ERPs:
% In fieldtrip epochs are called "trials"! In our information sampling
% (best-choice economic) task, every sequence sequence can have up to 10
% samples/epochs/trials(fieltrip).
% So, we define the conditions of the experimental task, and average ERPs
% as: 
% for every sequence, average 1:end-1 epochs and leave the last epoch as
% it is. This way, every sequence will have 2 ERPs (based on Nick's
% suggestions). This will give us 2 different cells of ERPs, one for
% all-but-last samples or for sampling-again and one for the last sample or
% price choice 


%% Define paths and variables 

basedir             = '/Users/christinadelta/Desktop/os_data';
task                = 'economic';
taskdir             = fullfile(basedir, task); 

subs                = dir(fullfile(taskdir, '*sub*'));
nsubs               = length(subs);

subname             = {subs.name}; % only keep the subject name 

% define initial variables 
blocks              = 2; 
blocktrials         = 20; 
totaltrials         = blocks * blocktrials;
nconds              = 2; % all-but-last (sample-again),last sample (choose price)

% define layout-montage for topoplots 
% check the biosemi64 layout that will be used for plotting topomaps
cfglayout           = [];
cfglayout.layout    = 'biosemi64.lay';
lo                  = ft_layoutplot(cfglayout);

%%  Load the data 

for subI = 1:nsubs
    
    % specify subI dir 
    subIdir = fullfile(taskdir, subname{subI});
    
    for blockI = 1:blocks
        
        % load subI block .bdf files 
        subFile = fullfile(subIdir, sprintf('sub_%02d_%s_block_%02d.bdf', subI, task, blockI));
        cfg                     = [];
        cfg.dataset             = subFile;

        cfg.trialdef.eventtype  = 'STATUS';
        cfg.trialdef.eventvalue = [1 100 101];
        cfg.trialdef.prestim    = 0.2;
        cfg.trialdef.poststim   = 0.8; 
        cfg                     = ft_definetrial(cfg);
        
       
        % this will help for timelock & tfr anlyses
        % extract trl list
        trl                     = cfg.trl;
        tstart                  = 100; % start of sequence
        tend                    = 101; % end of sequence
        
        % if the trl matrix starts with trigger code 101 - remove first
        % row
        if trl(1,4) == tend
            trl(1,:) = []; % remove entire first row 
        end

        blocksequences          = length(find(trl(:,4) == tstart));
        trialstart              = find(trl(:,4) == tstart);
        trialend                = find(trl(:,4) == tend);
               
        counter                 = 0;

        % loop over block-sequences 
        for iTrial = 1:blocksequences
            
            % add trial number 
            tmp                 = length(trialstart(iTrial)+1: trialend(iTrial)-1);
            for j = 1:tmp
                
                cnt             = counter + j; 
                trialnum(cnt,1) = iTrial;
                trialnum(cnt,2) = blockI;
            end 
            
            % update counter 
            counter             = counter + tmp ;
            
            clear tmp cnt
        end % end of sequence loop
        
        trl(trl(:,4) == tstart, :)  = [];
        trl(trl(:,4) == tend, :)    = [];
        trl(:,5)                    = trialnum(:,1); % add trialnum to the main list
        trl(:,6)                    = trialnum(:,2); % add trialnum to the main list
        
        clear trialend trialstart counter j trialnum
        
        % split the cfg data into all-but-last samples (sample-again) and last
        % sample (price choice)
        tmp_all             = 0;
        tmp_last            = 0;
        c                   = 0; % counter index
        l                   = 1; % last sample index
        
        for itrial = 1:blocksequences
            
            tmp                 = find(trl(:,5) == itrial);
            
            tl                  = length(tmp)-1; % trial length - 1 (these are the sample-again epochs)
            tmp_all(c+1:c+tl)   = tmp(1:end-1); % only pick the sample again epochs 
            tmp_last(:,l)       = tmp(end); % last epoch (price choice)
            
            % update c and l 
            c                   = c + tl;
            l                   = l + 1;
            
        end
        
        first_trl               = trl((tmp_all),:);
        last_trl                = trl((tmp_last),:);
        
        clear tmp_all tmp_last tmp itrial icondc l tl
        
        %% re-reference/preprocess

        cfg.reref           = 'yes';
        cfg.refchannel      = {'EXG1' 'EXG2'};
        cfg.demean          = 'yes';
        cfg.baselinewindow  = [-0.2 0];
        
        % re-write the trl list to the cfg struct and preprocess all data, easy/diff data & allButLast/last draws data structs 
        % the only thing that is changing here is "cfg.trl"
        cfg.trl             = trl;
        alldata             = ft_preprocessing(cfg);
        redata.alldata      = alldata; % update data struct
        clear cfg.trl
        
        cfg.trl             = first_trl;
        first_data          = ft_preprocessing(cfg);
        redata.first_data   = first_data;
        clear cfg.trl
        
        cfg.trl             = last_trl;
        last_data           = ft_preprocessing(cfg);
        redata.last_data    = last_data;
        
%         if redata.alldata.fsample ~= 512
%             cfg = [];
%             cfg.resamplefs = 512;
%             redata.alldata = ft_resampledata(cfg, redata.alldata);
%             
%             
%         end
%         
        % only keep eeg channels from now on % for all data structures
        % [exclude eog and mastoids]
        cfg                 = [];
        cfg.channel         = [1:64];
        
        data.alldata        = ft_selectdata(cfg, redata.alldata);
        data.firstdata      = ft_selectdata(cfg, redata.first_data);
        data.lastdata       = ft_selectdata(cfg, redata.last_data);
        
        % save preprocesssed block data in a .mat file 
        save(['economic_analysis/prepro/economic_preproc_sub_', num2str(subI), '_block_', num2str(blockI)], 'data')
        
        % clear data cfg 
        
        %% Optional: visualise the trials/epochs of interest 
        
        % COMMENT THIS SECTION IF YOU DON'T WANT TO VISUALISE THE EPOCHS 
        
        % NOTE: in typical analysis, each trial includes one stimulus and we epoch
        % the EEG signal around this stimulus. However, in the beads task each
        % trial contains a number of stimuli/draws (up to 10), thus here we define trials
        % and we epoch the EEG signal around each stimulus/draw in that sequnce.
        % Thus, we have 13 trials/sequences in each block, and length(trl) epochs. 

        % we can either visualise the epochs for each trial/sequence seperately or
        % all together. For simplicity, I'll visualise them together. 
%         cfg                 = [];
%         ft_databrowser(cfg, data)
        
    end % end of blocks loop
    
    % Load the matfiles, seperate them from the data struct and concatenate all blocks 
    % Load the matfiles and concatenate all blocks 
    for block = 1:blocks
        
        % if data structure is already in workspace, comment the part disp
        % and load parts 
        disp(['loading economic_analysis/prepro/economic_preproc_sub_', num2str(subI), '_block_', num2str(block)])
        load(['economic_analysis/prepro/economic_preproc_sub_', num2str(subI), '_block_', num2str(block)], 'data')
        
        partdata(block) = data.alldata;
        first_data(block)  = data.firstdata;
        last_data(block)   = data.lastdata;
        
        % visualise the block trials/epochs
%         cfg                 = [];
%         ft_databrowser(cfg, partdata(block))
        
    end
    
    % Append data
    cfg     = [];
    alldata = ft_appenddata(cfg, partdata(1), partdata(2));
    all_firstdata   = ft_appenddata(cfg, first_data(1), first_data(2));
    all_lastdata    = ft_appenddata(cfg, last_data(1), last_data(2));
    
   
    %% Remove artifacts with ICA 
    
    % THIS PART IS PROBABLY NOT NEEDED 
    % COMMENT IT IF ICA WONT RUN 
    % UNCOMMENT IF YOU WANT TO RUN ICA 
    
    % In general, in this experiment we don't need to run ICA for blink
    % removal because we instruct participants NOT to blink during this
    % time. Possibly will only need to filter.
    
    % Detect eog artifacts using ICA
%     cfg             = [];
%     cfg.method      = 'runica';
%     cfg.channel     = 1:64; % EEG channels only
%     datacomp        = ft_componentanalysis(cfg, alldata);
%     % save('beads_analysis/ica/datacomp', 'datacomp')
% 
%     % plot the components to detect the artifacts
%     figure
%     k = 1; f = 1;
% 
%     for icomp=1:length(datacomp.topo)
% 
%       if k>20
%         k = 1;
%         figure
%       end
% 
%       cfg           = [];
%       cfg.layout    = 'biosemi64.lay';
%       cfg.xlim      = [icomp icomp];
% 
%       subplot(4,5,k);
%       ft_topoplotER(cfg, datacomp);
%       title(icomp);
% 
%       k = k+1;
%     end
% 
%     % remove components that reflect eog artifacts
%     cfg             = [];
%     cfg.component   = [21 42 63 68]; % the exact numbers varies per run
%     clean_data      = ft_rejectcomponent(cfg, datacomp);
    
    %% filter data (clean) and run timelock (ERPs) analysis (4 analyses)
    
    % The ERPs that will mainly be used from now on are: [allButLastERPs, lastERPs]
    
    % Run ERPs analysis: 1 averaged ERP for sample 1:end-1 (first sample until last-1) & 1 averaged
    % ERP end (last sample- price choice) across all sequences/trials 
    % compute average ERPs for all-but-last samples
    cfg                             = [];
    cfg.preproc.lpfilter            = 'yes';
    cfg.preproc.lpfreq              = 40;
    cfg.trials                      = 'all';
    allButLastERPs                  = ft_timelockanalysis(cfg, all_firstdata);

    cfg                             = [];
    cfg.baseline                    = [-0.2 0];
    allButLastERPs                  = ft_timelockbaseline(cfg, allButLastERPs);
    
    % compute average ERPs for last draws
    cfg                             = [];
    cfg.preproc.lpfilter            = 'yes';
    cfg.preproc.lpfreq              = 40;
    cfg.trials                      = 'all';
    lastERPs                        = ft_timelockanalysis(cfg, all_lastdata);

    cfg                             = [];
    cfg.baseline                    = [-0.2 0];
    lastERPs                        = ft_timelockbaseline(cfg, lastERPs);
    
    % timelock analyses in cell for every subject
    allsubfirstERPs{1,subI}             = allButLastERPs;
    allsublastERPs{1,subI}              = lastERPs;
    
    % NOTE THAT IN THE TIMELOCK ANALYSIS ABOVE WE DID NOT INCLUDE THE
    % KEEPTRIALS OPTION WHICH MEANS THAT THE TRIALS ARE BEEING AVERAGED
    % (THIS WILL BE USED IN THE SECOND LEVEL/GROUP ANALYSIS) 
    % NOW WE CAN ALSO RUN TIMELOCK ANALYSIS AND SAVE ALL TRIALS INSTEAD OF
    % AVERAGED BY ADDING THE cfg.keeptrials = 'yes' OPTION - these will be
    % used for between trials stats (1st level)
    
    % Run ERPs analysis: 1 averaged ERP for samples 1:end-1 (first sample until last-1) & 1 averaged
    % ERP end (last sample/price choice) across all sequences/trials.
    % compute average ERPs for all-but-last samples
    cfg                             = [];
    cfg.preproc.lpfilter            = 'yes';
    cfg.preproc.lpfreq              = 40;
    cfg.keeptrials                  = 'yes';
    cfg.trials                      = 'all';
    allfirstERPs                    = ft_timelockanalysis(cfg, all_firstdata);

    cfg                             = [];
    cfg.baseline                    = [-0.2 0];
    allfirstERPs                    = ft_timelockbaseline(cfg, allfirstERPs);

    % compute average ERPs for last draws
    cfg                             = [];
    cfg.preproc.lpfilter            = 'yes';
    cfg.preproc.lpfreq              = 40;
    cfg.keeptrials                  = 'yes';
    cfg.trials                      = 'all';
    alllastERPs                     = ft_timelockanalysis(cfg, all_lastdata);

    cfg                             = [];
    cfg.baseline                    = [-0.2 0];
    alllastERPs                     = ft_timelockbaseline(cfg, alllastERPs);
    
    % save the ERPs in a cell for all subjects
    allsubs_allfirstERPs_keeptrials{1,subI} = allfirstERPs;
    allsubs_alllastERPs_keeptrials{1,subI}  = alllastERPs;
    
     %% Plot the ERPs
     
    % c) plot [all samples - last] vs last sample (over frontal
    % electrodes)
    cfg = [];
    cfg.channel = {'F1', 'F3', 'F5', 'F7', 'Fz', 'F2', 'F4', 'F6', 'F8'};
    figure; 
    ft_singleplotER(cfg, allButLastERPs, lastERPs)
    title('all draws but last vs last draw easy-frontal')
    legend('allButLast', 'last')
    
    % save figures
    print(gcf, '-dpng', ['economic_analysis/figures/erps/economic_sub', num2str(subI), '_frontal_allButLastERPs'])

    % d) plot [all draws - last] vs last draw for easy cond (over parietal
    % electrodes)
    cfg = [];
    cfg.channel = {'P1', 'P3', 'P5', 'P7', 'P9', 'Pz', 'P2', 'P4', 'P6', 'P8', 'P10'};
    figure; 
    ft_singleplotER(cfg, allButLastERPs, lastERPs)
    title('all draws but last vs last draw easy-parietal')
    legend('allButLast', 'last')
    
    print(gcf, '-dpng', ['economic_analysis/figures/erps/economic_sub', num2str(subI), '_parietal_allButLastERPs'])
  
    % ----------------------------------------------------------
    % plot ERPs on interactive mode (plot all but last and last)
    cfg             = [];
    cfg.layout      = 'biosemi64.lay';
    cfg.interactive = 'yes';
    cfg.showlabels  = 'yes';
    figure; ft_multiplotER(cfg, allButLastERPs, lastERPs)
    
    print(gcf, '-dpng', ['economic_analysis/figures/erps/economic_sub', num2str(subI), '_allchannels_ERPs'])

    % -----------------------------------------------------------------
    % maybe plot main ERPs in topoplots?
    % ALL but last
    cfg = [];
    cfg.xlim = [0.3 0.7];
    % cfg.zlim = [0 6e-14]; % this command messes up with the plot for some reason 
    cfg.layout  = 'biosemi64.lay';
    cfg.parameter = 'avg';
    cfg.interactive = 'yes';
    figure; ft_topoplotER(cfg, allButLastERPs); colorbar
    
    % save
    print(gcf, '-dpng', ['economic_analysis/figures/erps/economic_sub', num2str(subI), '_allbutlast_topo'])
    
    % Last, Easy
    cfg = [];
    cfg.xlim = [0.3 0.7];
    cfg.layout  = 'biosemi64.lay';
    cfg.parameter = 'avg';
    cfg.interactive = 'yes';
    figure; ft_topoplotER(cfg, lastERPs); colorbar
    
    print(gcf, '-dpng', ['economic_analysis/figures/erps/economic_sub', num2str(subI), '_last_topo'])
    
    %% Plot the differences (contrasts) using avg structures for each subject 
    
    % Plot the difference between all-but-last and last samples over frontal &
    % partietal sensors
    
    % first calculate the differences for all-but-last vs last samples ERPs
    cfg                     = [];
    cfg.operation           = 'x2-x1';
    cfg.parameter           = 'avg';
    difference_ERPs         = ft_math(cfg, lastERPs, allButLastERPs); % for easy trials 
   
    
    % plot over frontal sensors
    cfg = [];
    cfg.channel = {'F1', 'F3', 'F5', 'F7', 'Fz', 'F2', 'F4', 'F6', 'F8'};
    figure; 
    ft_singleplotER(cfg, difference_ERPs)
    title('forntal channels difference in all-but-last vs last ERPs')
    
    print(gcf, '-dpng', ['economic_analysis/figures/averages/economic_sub', num2str(subI), '_fig1_differenceERPallbutlast_vs_last'])
    
    % plot over parietal sensors 
    cfg = [];
    cfg.channel = {'P1', 'P3', 'P5', 'P7', 'P9', 'Pz', 'P2', 'P4', 'P6', 'P8', 'P10'};
    figure; 
    ft_singleplotER(cfg, difference_ERPs)
    title('parietal channels difference in all-but-last vs last ERPs')
    
    print(gcf, '-dpng', ['economic_analysis/figures/averages/economic_sub', num2str(subI), '_fig2_differenceERPallbutlast_vs_last'])
    
    % plot diferences as a movie
    % differences between all-but-last vs last sample ERPs 
    figure
    cfg        = [];
    cfg.layout = 'biosemi64.lay';
    ft_movieplotER(cfg, difference_ERPs); colorbar
    
    %% Time-Frequency Representation Analysis (TFR)
    
    % run TFR for Group 2 conditions (all-but-last samples vs last sample)
    cfg = [];
    cfg.channel    = 'all';
    cfg.method     = 'wavelet';
    cfg.width      = 7; % width definition is important 
    cfg.pad        = 10;
    cfg.output     = 'pow';
    cfg.foi        = 1:1:30;
    cfg.toi        = 'all'; % for computation efficiency use 'all'
    cfg.keeptrials = 'yes'; % keep the trials (will be used for statistical analysis)
    TFRwave_first  = ft_freqanalysis(cfg, all_firstdata);
    TFRwave_last   = ft_freqanalysis(cfg, all_lastdata);
    
    % save mat files in all subs cell 
    allsub_TFRfirst{1,subI}     = TFRwave_first;
    allsub_TFRlast{1,subI}      = TFRwave_last;
    
    %% Run TFR averaged for Within-Subjects analysis 
    
    % run TFR for Group 2 conditions (all-but-last draws vs last darw)
    cfg                 = [];
    cfg.channel         = 'all';
    cfg.method          = 'wavelet';
    cfg.width           = 7; % width definition is important 
    cfg.pad             = 10;
    cfg.output          = 'pow';
    cfg.foi             = 1:1:30;
    cfg.toi             = 'all'; % for computation efficiency use 'all'
    TFRwave_first_avg   = ft_freqanalysis(cfg, all_firstdata);
    TFRwave_last_avg    = ft_freqanalysis(cfg, all_lastdata);
    
    % save mat files in all subs cell 
    allsub_TFRfirst_avg{1,subI}     = TFRwave_first_avg;
    allsub_TFRlast_avg{1,subI}      = TFRwave_last_avg;
    
    cfg = [];
    cfg.baseline     = [-0.2 0];
    cfg.baselinetype = 'absolute';
    cfg.marker       = 'on';
    cfg.showlabels   = 'yes';
    cfg.layout       = 'biosemi64.lay';
    figure; ft_multiplotTFR(cfg, TFRwave_first_avg); colorbar; ft_multiplotTFR(cfg, TFRwave_last_avg); colorbar;
    
    % single plots 
    cfg = [];
    cfg.baseline     = [-0.2 0];
    cfg.baselinetype = 'absolute';
    cfg.marker       = 'on';
    cfg.maskstyle    = 'saturation';
    cfg.channel      = {'P1', 'P3', 'P5', 'P7', 'P9', 'Pz', 'P2', 'P4', 'P6', 'P8', 'P10'};
    cfg.layout       = 'biosemi64.lay';
    
    figure; ft_singleplotTFR(cfg, TFRwave_first_avg); colorbar;
    print(gcf, '-dpng', ['economic_analysis/figures/tfr/economic_sub', num2str(subI), '_tfr_all_but_last'])
    
    figure; ft_singleplotTFR(cfg, TFRwave_last_avg); colorbar;
    print(gcf, '-dpng', ['economic_analysis/figures/tfr/economic_sub', num2str(subI), '_tfr_last'])
    

end % end of subjects loop

%% Save all sub cell files for later trial-level & group-level analyses
    
% save all subjects ERP analyes 
save('economic_analysis/erps/economic_allsubsAllButLastERPs', 'allsubfirstERPs')
save('economic_analysis/erps/economic_allsubsLastERPs', 'allsublastERPs')

save('economic_analysis/erps/economic_allsubs_keeptrialsFirstERPs', 'allsubs_allfirstERPs_keeptrials')
save('economic_analysis/erps/economic_allsubs_keeptrialsLastERPs', 'allsubs_alllastERPs_keeptrials')

% save all subjects TFR analyes - with the full trials
save('economic_analysis/tfr/economic_allsubs_tfrFirst', 'allsub_TFRfirst', '-v7.3')
save('economic_analysis/tfr/economics_allsubs_tfrLast', 'allsub_TFRlast', '-v7.3')

% save all subject TFR analyses - averages
save('economic_analysis/tfr/economic_allsubs_tfrFirst_avg', 'allsub_TFRfirst_avg', '-v7.3')
save('economic_analysis/tfr/economic_allsubs_tfrLast_avg', 'allsub_TFRlast_avg', '-v7.3')

%% Compute TFR grand averages for all datasets/conditions

% first load TFRs if not loaded 
% load all subjects TFR analyes - with the full trials
load('economic_analysis/tfr/economic_allsubs_tfrFirst', 'allsub_TFRfirst')
load('economic_analysis/tfr/economics_allsubs_tfrLast', 'allsub_TFRlast')

% load all subject TFR analyses - averages
load('economic_analysis/tfr/economic_allsubs_tfrFirst_avg', 'allsub_TFRfirst_avg')
load('economic_analysis/tfr/economic_allsubs_tfrLast_avg', 'allsub_TFRlast_avg')

% calculate grand average for each condition with the individual data set
% to yes
cfg                 = [];
cfg.channel         = 'all';
cfg.toilim          = 'all';
cfg.foilim          = 'all';
cfg.keepindividual  = 'yes';
cfg.parameter       = 'powspctrm';

ga_TFRfirst_ind     = ft_freqgrandaverage(cfg, allsub_TFRfirst_avg{:});
ga_TFRlast_ind      = ft_freqgrandaverage(cfg, allsub_TFRlast_avg{:});

% calculate grand average for each condition 
cfg                 = [];
cfg.channel         = 'all';
cfg.toilim          = 'all';
cfg.foilim          = 'all';
cfg.parameter       = 'powspctrm';

ga_TFRfirst         = ft_freqgrandaverage(cfg, allsub_TFRfirst_avg{:});
ga_TFRlast          = ft_freqgrandaverage(cfg, allsub_TFRlast_avg{:});

% save TFR grand averages (with individual trials) 
save('economic_analysis/tfr/economic_ga_TFRfirst_ind', 'ga_TFRfirst_ind', '-v7.3')
save('economic_analysis/tfr/economic_ga_TFRlast_ind', 'ga_TFRlast_ind', '-v7.3')

% save TFR grand averages (without individual trials) 
save('economic_analysis/tfr/economic_ga_TFRfirst', 'ga_TFRfirst', '-v7.3')
save('economic_analysis/tfr/economic_ga_TFRlast', 'ga_TFRlast', '-v7.3')

%% Compute Timelock Grand Averages for each dataset/condition 

% load all subjects ERP analyes if needed 
load('economic_analysis/erps/economic_allsubsAllButLastERPs', 'allsubfirstERPs')
load('economic_analysis/erps/economic_allsubsLastERPs', 'allsublastERPs')

load('economic_analysis/erps/economic_allsubs_keeptrialsFirstERPs', 'allsubs_allfirstERPs_keeptrials')
load('economic_analysis/erps/economic_allsubs_keeptrialsLastERPs', 'allsubs_alllastERPs_keeptrials')

% now calculate timelocked grand averages 
cfg             = [];
cfg.channel     = 'all';
cfg.latency     = 'all';
cfg.parameter   = 'avg';

gaERP_first     = ft_timelockgrandaverage(cfg, allsubfirstERPs{:});
gaERP_last      = ft_timelockgrandaverage(cfg, allsublastERPs{:});

% save timelock averages 
save('economic_analysis/erps/economic_gaERP_first', 'gaERP_first')
save('economic_analysis/erps/economic_gaERP_last', 'gaERP_last')

% load the ERP grand averages if needed 
% load('economic_analysis/erps/economic_gaERP_first', 'gaERP_first')
% load('economic_analysis/erps/economic_gaERP_last', 'gaERP_last')

%% Visualise grand average ERPs and compute GA contrasts

% plot grand averages for first vs last ERPs over frontal sites 
cfg = [];
cfg.channel = {'F1', 'F3', 'F5', 'F7', 'Fz', 'F2', 'F4', 'F6', 'F8'};
figure; 
ft_singleplotER(cfg, gaERP_first, gaERP_last)
title('forntal channels GA ERPs')
legend('first', 'last')

print(gcf, '-dpng', 'economic_analysis/figures/grand_averageERPs_plots/economic_GA_first_vs_last_ERP_frontal')

% ------------------------------

% plot grand averages for first vs last ERPs over parietal sites 
cfg = [];
cfg.channel = {'P1', 'P3', 'P5', 'P7', 'P9', 'Pz', 'P2', 'P4', 'P6', 'P8', 'P10'};
figure; 
ft_singleplotER(cfg, gaERP_first, gaERP_last)
title('Parietal channels GA ERPs')
legend('first', 'last')

print(gcf, '-dpng', 'economic_analysis/figures/grand_averageERPs_plots/economic_GA_first_vs_last_ERP_parietal')

%% Compute GA differences for each codnition/dataset

% first compute the averages 
cfg                         = [];
cfg.operation               = 'x2-x1';
cfg.parameter               = 'avg';
difference_lastvsfirst      = ft_math(cfg, gaERP_last, gaERP_first);

% plot difference between last vs first over frontal sites 
cfg = [];
cfg.channel = {'F1', 'F3', 'F5', 'F7', 'Fz', 'F2', 'F4', 'F6', 'F8'};
figure; 
ft_singleplotER(cfg, difference_lastvsfirst)
title('GA difference first vs last ERPs over frontal sites')

print(gcf, '-dpng', 'economic_analysis/figures/differences/economic_figGA_difference_lastvsfirst_frontal')

% plot difference between easy vs difficult over frontal sites 
cfg = [];
cfg.channel = {'P1', 'P3', 'P5', 'P7', 'P9', 'Pz', 'P2', 'P4', 'P6', 'P8', 'P10'};
figure; 
ft_singleplotER(cfg, difference_lastvsfirst)
title('GA difference first vs last ERPs over parietal sites')

print(gcf, '-dpng', 'economic_analysis/figures/differences/economic_figGA_difference_lastvsfirst_parietal')


%% Visualise grand average TFRs for all conditions/datasets

%load GA TFRs if needed 
% load('economic_analysis/tfr/economic_ga_TFRfirst_ind', 'ga_TFRfirst_ind')
% load('economic_analysis/tfr/economic_ga_TFRlast_ind', 'ga_TFRlast_ind')
% load('economic_analysis/tfr/economic_ga_TFRfirst', 'ga_TFRfirst')
% load('economic_analysis/tfr/economic_ga_TFRlast', 'ga_TFRlast')

% plot GA TFR for first vs last epochs over all sensors 
cfg                 = [];
cfg.baseline        = [-0.2 0];
cfg.baselinetype    = 'absolute';
cfg.marker          = 'on';
cfg.showlabels      = 'yes';
cfg.layout          = 'biosemi64.lay';
figure; ft_multiplotTFR(cfg, ga_TFRfirst); colorbar; ft_multiplotTFR(cfg, ga_TFRlast); colorbar;

% single plot of TFR first vs last over parietal sites
cfg = [];
cfg.baseline     = [-0.2 0];
cfg.baselinetype = 'absolute';
cfg.marker       = 'on';
cfg.maskstyle    = 'saturation';
cfg.channel      = {'P1', 'P3', 'P5', 'P7', 'P9', 'Pz', 'P2', 'P4', 'P6', 'P8', 'P10'};
cfg.layout       = 'biosemi64.lay';

figure; ft_singleplotTFR(cfg, ga_TFRfirst); colorbar;
print(gcf, '-dpng', 'economic_analysis/figures/grand_averageTFRs_plots/economic_tfr_GA_first')

figure; ft_singleplotTFR(cfg, ga_TFRlast); colorbar;
print(gcf, '-dpng', 'economic_analysis/figures/grand_averageTFRs_plots/economic_tfr_GA_last')

%% Run ERP statistics (between-trials/first level)

% load all subs cells if needed
load('economic_analysis/erps/economic_allsubs_keeptrialsFirstERPs', 'allsubs_allfirstERPs_keeptrials')
load('economic_analysis/erps/economic_allsubs_keeptrialsLastERPs', 'allsubs_alllastERPs_keeptrials')

% run between trial statistics  

%% Within-Subjects statitsics ERPs 

% load all subjects ERP analyes if needed 
load('economic_analysis/erps/economic_allsubsAllButLastERPs', 'allsubfirstERPs')
load('economic_analysis/erps/economic_allsubsLastERPs', 'allsublastERPs')

%% Run group stats to first vs last trials 

% use the original data of one sub as example dataset to feed the cfg.neighbours matrix
% (to get sensor labels)
tmp_sub     = 2;
blocks      = 2;


for block = 1:blocks
    
    % if data structure is already in workspace, comment the part disp
    % and load parts 
    disp(['loading economic_analysis/prepro/economic_preproc_sub_', num2str(tmp_sub), '_block_', num2str(block)])
    load(['economic_analysis/prepro/economic_preproc_sub_', num2str(tmp_sub), '_block_', num2str(block)], 'data')

    first_data(block)  = data.firstdata;
    last_data(block)   = data.lastdata;

end

% Append data
cfg             = [];
all_firstdata   = ft_appenddata(cfg, first_data(1), first_data(2));
all_lastdata    = ft_appenddata(cfg, last_data(1), last_data(2));

subs            = length(allsubfirstERPs);

% first create the neighbours structure 
cfg_neighb              = [];
cfg_neighb.layout       = 'biosemi64.lay'; %in meters
cfg_neighb.method       = 'distance';
neighbours              = ft_prepare_neighbours(cfg_neighb, all_firstdata);

cfg                     = [];
cfg.neighbours          = neighbours;   

cfg.channel             = 'all';        
cfg.latency             = [0.3 0.8];   


% permutation tests
cfg.method              = 'montecarlo';    
cfg.statistic           = 'depsamplesT'; 

cfg.correctm            = 'cluster';
cfg.clusteralpha        = 0.05;         

cfg.clusterstatistic    = 'maxsum';                            
cfg.minnbchan           = 2;            
cfg.tail                = 0;           
cfg.clustertail         = 0;
cfg.alpha               = 0.025;        
cfg.numrandomization    = 1000;       

% create the design matrix
design                  = zeros(2, subs*2);
design(1,:)             = [1:subs 1:subs];
design(2,:)             = [ones(1,subs) ones(1,subs)*2];

cfg.design              = design;
cfg.uvar                = 1; % row 1 of design matrix contain var 1
cfg.ivar                = 2; % row 2 of design matrix contain var 2

[stat] = ft_timelockstatistics(cfg, allsublastERPs{:}, allsubfirstERPs{:});

% save stats 
save('economic_analysis/erps_group_stats/economic_stats_group_last_vs_first', 'stat')

% Make a vector of all p-values associated with the clusters from ft_timelockstatistics.
pos_cluster_pvals   = [stat.posclusters(:).prob];

% Then, find which clusters are deemed interesting to visualize, here we use a cutoff criterion based on the
% cluster-associated p-value, and take a 5% two-sided cutoff (i.e. 0.025 for the positive and negative clusters,
% respectively
pos_clust           = find(pos_cluster_pvals < 0.025);
pos                 = ismember(stat.posclusterslabelmat, pos_clust);


% and now for the negative clusters...
neg_cluster_pvals   = [stat.negclusters(:).prob];
neg_clust           = find(neg_cluster_pvals < 0.025);
neg                 = ismember(stat.negclusterslabelmat, neg_clust);

% if any positive clusters survived the cutoff.. loop over all sig positive clusters
if ~isempty(pos_clust)
    for i=pos_clust

        cfg=[];
        cfg.highlight = 'on';
        cfg.zparam    = 'stat';
        cfg.layout    = 'biosemi64.lay';
        cfg.style     = 'straight';
        cfg.gridscale = 500;

        % find the significant time range for this cluster
        tmp=[];
        for t = 1:length(stat.time)
            if ~isempty(find(any(stat.posclusterslabelmat(:,t)==pos_clust)))
              tmp = [tmp t];
            end
        end
        cfg.xlim      = [stat.time(tmp(1)) stat.time(tmp(end))];

        % find the channels belonging to this cluster
        cfg.highlightchannel = [];

        for c = 1:length(stat.label)
            if ~isempty(find(any(stat.posclusterslabelmat(:, c)==pos_clust)))
              cfg.highlightchannel = [cfg.highlightchannel c];
            end
        end

        figure
        ft_topoplotER(cfg, stat);
        title('positive cluster')
        print(gcf, '-dpng', ['economic_analysis/figures/group_stats_erps/fig1_STAT_pos', num2str(i)])
    end
end

% now run the same for negative clusters 
if ~isempty(neg_clust)

    for i=neg_clust

        cfg=[];
        cfg.highlight = 'on';
        cfg.zparam    = 'stat';
        cfg.layout    = 'biosemi64.lay';
        cfg.style     = 'straight';
        cfg.gridscale = 500;

        % find the significant time range for this cluster
        tmp=[];

        for t = 1:length(stat.time)
            if ~isempty(find(any(stat.negclusterslabelmat(:,t)==neg_clust)))
              tmp = [tmp t];
            end
        end
        cfg.xlim      = [stat.time(tmp(1)) stat.time(tmp(end))];

        % find the channels belonging to this cluster
        cfg.highlightchannel = [];
        for c = 1:length(stat.label)
            if ~isempty(find(any(stat.negclusterslabelmat(:,c)==neg_clust)))
              cfg.highlightchannel = [cfg.highlightchannel c];
            end
        end

        figure
        ft_topoplotER(cfg, stat);
        title('negative cluster')
        print(gcf, '-dpng', ['economic_analysis/figures/group_stats_erps/fig2_STAT_neg', num2str(i)])
    end
end

% try ploting on clusterplots for time-windows, more informative than the plots above!!
% negative & positive clusters
% find max and min t-values first 
if ~isempty(pos_clust) | ~isempty(neg_clust) 
    maxval                  = max(stat.stat, [], 'all');
    minval                  = min(stat.stat, [], 'all');

    cfg                     = [];
    cfg.zlim                = [minval maxval]; % T-values
    cfg.alpha               = 0.025;
    cfg.highlightcolorpos   = [0 0 0];
    cfg.highlightcolorneg   = [1 1 0];
    cfg.saveaspng           = 'economic_analysis/figures/group_stats_erps/economic_group_last_vs_first';
    cfg.toi                 = 0.3:0.01:0.6; % times of interest
    cfg.layout              = 'biosemi64.lay';
    ft_clusterplot(cfg,stat); 
end

clear stat cfg cfg_neighb neigbours design

%% Run TFR statistics (between-trials/first level)


% load stuff
% run stats



%% Run TFR statistics (Within-Subjects/second level)

%load GA TFRs if needed 
% load('economic_analysis/tfr/economic_ga_TFRfirst_ind', 'ga_TFRfirst_ind')
% load('economic_analysis/tfr/economic_ga_TFRlast_ind', 'ga_TFRlast_ind')
% load('economic_analysis/tfr/economic_ga_TFRfirst', 'ga_TFRfirst')
% load('economic_analysis/tfr/economic_ga_TFRlast', 'ga_TFRlast')

%% Run TFR statistics for last vs first trials

% create the neighbours structure which will be critical to compute stats
cfg_neighb              = [];
cfg_neighb.layout       = 'biosemi64.lay'; %in meters
cfg_neighb.method       = 'distance';
neighbours              = ft_prepare_neighbours(cfg_neighb, ga_TFRfirst_ind);

% permutations configuration structure 
cfg                     = [];
cfg.neighbours          = neighbours;   % the neighbours specify for each sensor with
                                        % which other sensors it can form clusters
cfg.channel             = 'all';        % cell-array with selected channel labels
cfg.latency             = [0 0.8];      % time interval over which the experimental - I will use 'all' for now
                                        % conditions must be compared (in seconds)

cfg.frequency           = 30;           % frequences of interest 

% permutation tests
cfg.method              = 'montecarlo';     % use the Monte Carlo Method to calculate the significance probability
cfg.statistic           = 'ft_statfun_depsamplesT';  % use the dependent samples T-statistic as a measure to
                                            % evaluate the effect at the group level
cfg.correctm            = 'cluster';
cfg.clusteralpha        = 0.05;         % alpha level of the sample-specific test statistic that
                                        % will be used for thresholding
cfg.clusterstatistic    = 'maxsum';     % test statistic that will be evaluated under the
                                        % permutation distribution.
cfg.minnbchan           = 2;            % minimum number of neighborhood channels that is
                                        % required for a selected sample to be included
                                        % in the clustering algorithm (default=0).
cfg.tail                = 0;            % -1, 1 or 0 (default = 0); one-sided or two-sided test
cfg.clustertail         = 0;
cfg.alpha               = 0.025;        % alpha level of the permutation test (if two-sided set to 0.025)
cfg.numrandomization    = 1000;         % number of draws from the permutation distribution

% WITHIN SUBJECTS DESIGN MATRIX
subj                    = 18; % this needs to be changed everytime we add a subject
design                  = zeros(2,2*subj);

for i = 1:subj
  design(1,i)           = i;
end

for i = 1:subj
  design(1,subj+i)      = i;
end

design(2,1:subj)        = 1;
design(2,subj+1:2*subj) = 2;

cfg.design              = design;
cfg.uvar                = 1;
cfg.ivar                = 2;

[stat] = ft_freqstatistics(cfg, ga_TFRlast_ind, ga_TFRfirst_ind);

% save stats 
save('economic_analysis/tfr_group_stats/economic_TFRstats_group_last_vs_first', 'stat')

% Make a vector of all p-values associated with the clusters from ft_timelockstatistics.
pos_cluster_pvals   = [stat.posclusters(:).prob];

% Then, find which clusters are deemed interesting to visualize, here we use a cutoff criterion based on the
% cluster-associated p-value, and take a 5% two-sided cutoff (i.e. 0.025 for the positive and negative clusters,
% respectively
pos_clust           = find(pos_cluster_pvals < 0.025);
pos                 = ismember(stat.posclusterslabelmat, pos_clust);


% and now for the negative clusters...
neg_cluster_pvals   = [stat.negclusters(:).prob];
neg_clust           = find(neg_cluster_pvals < 0.025);
neg                 = ismember(stat.negclusterslabelmat, neg_clust);

% if there is a cluster that survived the cutoff, plot the result 
% cfg = [];
% cfg.alpha  = 0.025; % if alpha = 0.025 doesn't work try 0.05
% cfg.parameter = 'stat';
% cfg.zlim   = [-4 4];
% cfg.layout = 'biosemi64.lay';
% ft_clusterplot(cfg, stat);

% try ploting on clusterplots for time-windows, more informative than the plots above!!
% negative & positive clusters
% find max and min t-values first 
if ~isempty(pos_clust) | ~isempty(neg_clust) 
    maxval                  = max(stat.stat, [], 'all');
    minval                  = min(stat.stat, [], 'all');

    cfg                     = [];
    cfg.zlim                = [minval maxval]; % T-values
    cfg.alpha               = 0.025;
    cfg.highlightcolorpos   = [0 0 0];
    cfg.highlightcolorneg   = [1 1 0];
    cfg.saveaspng           = 'economic_analysis/figures/group_stats_tfrs/economic_groupTFR_last_vs_first';
    %cfg.toi                 = 0.0:0.01:0.6; % times of interest
    cfg.layout              = 'biosemi64.lay';
    ft_clusterplot(cfg,stat); 
end

clear stat cfg cfg_neighb






