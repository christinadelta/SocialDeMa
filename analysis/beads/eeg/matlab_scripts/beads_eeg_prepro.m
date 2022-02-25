%% BEADS FORMAL PREPROCESSING AND ANALYSES SCRIPT W/ FIELDTRIP

% christinadelta 
% Date: Jan 2022
% last Update 20/02/2022

% final preprocessing and analyses of the EEG data - beads task using
% FIELDTRIP 

% A few differences in pre-processing methods between biosemi (.bdf) data
% and other devises is that:
% 1) with biosemi with don't need to do re-referencing to obtain EOG channels,
% the signal is recorded as bipolar 
% 2) biosemi records and saves the signal UN-REFERENCED, thus, we don't
% include an "implicitref" during analysis. We only need to provide the "reref" channel(s)

% HOW DO WE COMPUTE ERPs:
% In fieldtrip epochs are called "trials"! In our information sampling
% (beads) task, every sequence sequence can have up to 10
% samples/epochs/trials(fieltrip).
% So, we define the conditions of the experimental task, and average ERPs
% at 3 different levels:
% 1) average ERPs for each condition across sequences 
% 2) don't average ERPs (just epoch the signal based on the conditions) and
% leave ERPs as they are (without averaging
% 3) for every sequence, average 1:end-1 epochs and leave the last epoch as
% it is. This way, every sequence will have 2 ERPs (based on Nick's
% suggestions)

% TODO:
% 1). Visualise ERPs in different ways (done)
% 2). Run stats
% 3). Extract Frequences (done)

%% Define paths and variables 

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
nconds              = 4; % easy blue, easy green, diff blue, diff green
totalconds          = 2; % easy, diff
conditions          = {'Easy', 'Difficult'};

% define layout-montage for topoplots 
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
        
        % load subI block .bdf files 
        subFile = fullfile(subIdir, sprintf('sub_%02d_%s_block_%02d.bdf', subI, task, blockI));
        cfg                     = [];
        cfg.dataset             = subFile;

        cfg.trialdef.eventtype  = 'STATUS';
        cfg.trialdef.eventvalue = [1 2 3 4 102 103];
        cfg.trialdef.prestim    = 0.2;
        cfg.trialdef.poststim   = 0.8; 
        cfg                     = ft_definetrial(cfg);
        
        %% Split block-trials in different preprocessing data structures 
        
        % this will help for timelock & tfr anlyses
        % extract trl list
        trl                     = cfg.trl;
        tstart                  = 102; % start of sequence
        tend                    = 103; % end of sequence

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

        % remove trialstart and trialend from the trl list
        trl(trl(:,4) == tstart, :)  = [];
        trl(trl(:,4) == tend, :)    = [];
        trl(:,5)                    = trialnum(:,1); % add trialnum to the main list
        trl(:,6)                    = trialnum(:,2); % add trialnum to the main list
        
        % now split data to conditions for later preprocessing seperately, & move condition data to
        % new data structures     
        trl_length                      = length(trl);

        for i = 1:trl_length

            if trl(i,4) == 1 | trl(i,4) == 2
                trl(i,7) = 1;

            elseif trl(i,4) == 3 | trl(i,4) == 4
                trl(i,7) = 2;

            end
        end
        
        % clear for memory eficiency
        clear trialend trialstart j i counter cnt trialnum 

        % split cfg data structure in easy and difficult conditions 
        trl_easy            = trl((find(trl(:,7) == 1)),:);
        trl_diff            = trl((find(trl(:,7) == 2)),:);
        
        % split cfg data struct into [all draws - last] and [last draw]
        tmp_all             = 0;
        tmp_last            = 0;
        c                   = 0; % counter index
        l                   = 1; % last draw index
        
        for icond = 1:totalconds
            for itrial = 1:blocktrials

                tmp = find(trl(:,7)== icond & trl(:,5)== itrial);

                if ~isempty(tmp)
                    tl                  = length(tmp)-1;
                    tmp_all(c+1:c+tl)   = tmp(1:end-1); % only pick 
                    tmp_last(:,l)       = tmp(end);

                    % update c and l 
                    c                   = c + tl;
                    l                   = l + 1;
                end % end of if
            end % end of trials loop
        end % end of iconds loop
        
        first_trl                      = trl((tmp_all),:);
        last_trl                       = trl((tmp_last),:);
        
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
        
        cfg.trl             = trl_easy;
        easy_data           = ft_preprocessing(cfg);
        redata.easy_data    = easy_data;
        clear cfg.trl
        
        cfg.trl             = trl_diff;
        diff_data           = ft_preprocessing(cfg);
        redata.diff_data    = diff_data;
        clear cfg.trl
        
        cfg.trl             = first_trl;
        first_data          = ft_preprocessing(cfg);
        redata.first_data   = first_data;
        clear cfg.trl
        
        cfg.trl             = last_trl;
        last_data           = ft_preprocessing(cfg);
        redata.last_data    = last_data;
        
        % only keep eeg channels from now on % for all data structures
        % [exclude eog and mastoids]
        cfg                 = [];
        cfg.channel         = [1:64];
        
        data.alldata        = ft_selectdata(cfg, redata.alldata);
        data.easydata       = ft_selectdata(cfg, redata.easy_data);
        data.diffdata       = ft_selectdata(cfg, redata.diff_data);
        data.firstdata      = ft_selectdata(cfg, redata.first_data);
        data.lastdata       = ft_selectdata(cfg, redata.last_data);
        
        % save preprocesssed block data in a .mat file 
        save(['beads_analysis/prepro/beads_preproc_sub_', num2str(subI), '_block_', num2str(blockI)], 'data')
        
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
    
    %% Append data into one structure

    % clear all
    
    % Load the matfiles, seperate them from the data struct and concatenate all blocks 
    % Load the matfiles and concatenate all blocks 
    for block = 1:blocks
        
        % if data structure is already in workspace, comment the part disp
        % and load parts 
        disp(['loading beads_analysis/prepro/beads_preproc_sub_', num2str(subI), '_block_', num2str(block)])
        load(['beads_analysis/prepro/beads_preproc_sub_', num2str(subI), '_block_', num2str(block)], 'data')
        
        partdata(block) = data.alldata;
        first_data(block)  = data.firstdata;
        last_data(block)   = data.lastdata;
        
        % visualise the block trials/epochs
%         cfg                 = [];
%         ft_databrowser(cfg, partdata(block))
        
    end
    
    % Append data
    cfg     = [];
    alldata = ft_appenddata(cfg, partdata(1), partdata(2), partdata(3), partdata(4));
    all_firstdata   = ft_appenddata(cfg, first_data(1), first_data(2), first_data(3), first_data(4));
    all_lastdata    = ft_appenddata(cfg, last_data(1), last_data(2), last_data(3), last_data(4));
    
    
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
%     datacomp        = ft_componentanalysis(cfg, data);
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
    
    % I will run 3 different ERP analyses only to be able to visualise them
    % properly and to explore the data. 
    
    % The ERPs that will mainly be used from now on are: [avgERPtwocond, allButLastERPs, lastERPs]
    
    % 1) Run ERPs analysis: 1 average ERP for each sequence
    % Average ERP: 1:end (epochs) ( 52 averaged ERPs total) - 4 Conditions
    for icond = 1:nconds 
        cfg                             = [];
        cfg.preproc.lpfilter            = 'yes';
        cfg.preproc.lpfreq              = 40;
        cfg.trials                      = find(alldata.trialinfo(:,1) == icond);
        avgERPallcond{icond}            = ft_timelockanalysis(cfg, alldata);
        
        cfg                             = [];
        cfg.baseline                    = [-0.2 0];
        avgERPallcond{icond}            = ft_timelockbaseline(cfg, avgERPallcond{icond});
    end
    
    % 2) Run ERPs analysis: 1 average ERP for each sequence
    % Average ERP: 1:end (epochs) ( 52 averaged ERPs total) - 2 Conditions
    
    % add a few more variables 
    totalconds      = 2; 
    totaldraws      = length(alldata.trialinfo);
    condtrials      = blocktrials * totalconds;
    
    % run 2ND ERP analysis using 2 conditions (easy vs diff for each
    % sequence)
    for icond = 1:totalconds 
        
        cfg                             = [];
        cfg.preproc.lpfilter            = 'yes';
        cfg.preproc.lpfreq              = 40;
        cfg.trials                      = find(alldata.trialinfo(:,4) == icond);
        avgERPtwocond{icond}            = ft_timelockanalysis(cfg, alldata);
        
        cfg                             = [];
        cfg.baseline                    = [-0.2 0];
        avgERPtwocond{icond}            = ft_timelockbaseline(cfg, avgERPtwocond{icond});
        
        
    end   

    % 3) Run ERPs analysis: 1 averaged ERP for draws 1:end-1 (first draw until last-1) & 1 averaged
    % ERP end (last draw) across all sequences/trials for each condition
    % (easy, difficult). For this we first need to make
    % sure that for each condition (easy, diff) we get 2 ERPs; so, let's
    % first split the draws in different matrices. One that will contain
    % all the [1: end-1] draws and one matrix that will contain all the
    % [end] draws. Then run timelock analysis (ERPs) and store in cell for
    % each condition separately.
    
    for icond = 1:totalconds
        
        % compute average ERPs for all-but-last draws
        cfg                             = [];
        cfg.preproc.lpfilter            = 'yes';
        cfg.preproc.lpfreq              = 40;
        cfg.trials                      = find(all_firstdata.trialinfo(:,4) == icond);
        allButLastERPs{icond}           = ft_timelockanalysis(cfg, all_firstdata);
        
        cfg                             = [];
        cfg.baseline                    = [-0.2 0];
        allButLastERPs{icond}           = ft_timelockbaseline(cfg, allButLastERPs{icond});
        
        % compute average ERPs for last draws
        cfg                             = [];
        cfg.preproc.lpfilter            = 'yes';
        cfg.preproc.lpfreq              = 40;
        cfg.trials                      = find(all_lastdata.trialinfo(:,4) == icond);
        lastERPs{icond}                 = ft_timelockanalysis(cfg, all_lastdata);
        
        cfg                             = [];
        cfg.baseline                    = [-0.2 0];
        lastERPs{icond}                 = ft_timelockbaseline(cfg, lastERPs{icond});
         
    end
    
    % save timelock analyses in cell for every subject
    allsubERPallcond{1,subI}            = avgERPallcond;
    allsubERPtwocond{1,subI}            = avgERPtwocond;
    allsubfirstERPs{1,subI}             = allButLastERPs;
    allsublastERPs{1,subI}              = lastERPs;
   
    % NOTE THAT IN THE TIMELOCK ANALYSIS ABOVE WE DID NOT INCLUDE THE
    % KEEPTRIALS OPTION WHICH MEANS THAT THE TRIALS ARE BEEING AVERAGED
    % (THIS IS HELPFUL FOR COMPUTING DIFFERENCE BETWEEN THE CONDITIONS) 
    % NOW WE CAN ALSO RUN TIMELOCK ANALYSIS AND SAVE ALL TRIALS INSTEAD OF
    % AVERAGED BY ADDING THE cfg.keeptrials = 'yes' OPTION - these will be
    % used for between trials stats (1st level)
    
    % 1) Run ERPs analysis: 4 Conditions
    for icond = 1:nconds 
        cfg                             = [];
        cfg.preproc.lpfilter            = 'yes';
        cfg.preproc.lpfreq              = 40;
        cfg.keeptrials                  = 'yes';
        cfg.trials                      = find(alldata.trialinfo(:,1) == icond);
        ERPallcond{icond}               = ft_timelockanalysis(cfg, alldata);
        
        cfg                             = [];
        cfg.baseline                    = [-0.2 0];
        ERPallcond{icond}               = ft_timelockbaseline(cfg, ERPallcond{icond});
    end
    
    % 2) Run ERPs analysis: 2 Conditions
    for icond = 1:totalconds 
        
        cfg                             = [];
        cfg.preproc.lpfilter            = 'yes';
        cfg.preproc.lpfreq              = 40;
        cfg.keeptrials                  = 'yes';
        cfg.trials                      = find(alldata.trialinfo(:,4) == icond);
        ERPtwocond{icond}               = ft_timelockanalysis(cfg, alldata);
        
        cfg                             = [];
        cfg.baseline                    = [-0.2 0];
        ERPtwocond{icond}               = ft_timelockbaseline(cfg, ERPtwocond{icond});
        
        
    end   

    % 3) Run ERPs analysis: 1 averaged ERP for draws 1:end-1 (first draw until last-1) & 1 averaged
    % ERP end (last draw) across all sequences/trials for each condition
    % (easy, difficult).
    for icond = 1:totalconds
        
        % compute average ERPs for all-but-last draws
        cfg                             = [];
        cfg.preproc.lpfilter            = 'yes';
        cfg.preproc.lpfreq              = 40;
        cfg.keeptrials                  = 'yes';
        cfg.trials                      = find(all_firstdata.trialinfo(:,4) == icond);
        allfirstERPs{icond}             = ft_timelockanalysis(cfg, all_firstdata);
        
        cfg                             = [];
        cfg.baseline                    = [-0.2 0];
        allfirstERPs{icond}             = ft_timelockbaseline(cfg, allfirstERPs{icond});
        
        % compute average ERPs for last draws
        cfg                             = [];
        cfg.preproc.lpfilter            = 'yes';
        cfg.preproc.lpfreq              = 40;
        cfg.keeptrials                  = 'yes';
        cfg.trials                      = find(all_lastdata.trialinfo(:,4) == icond);
        alllastERPs{icond}              = ft_timelockanalysis(cfg, all_lastdata);
        
        cfg                             = [];
        cfg.baseline                    = [-0.2 0];
        alllastERPs{icond}              = ft_timelockbaseline(cfg, alllastERPs{icond});
         
    end
    
    % save the ERPs in a cell for all subjects
    allsubs_ERPallcond_keeptrials{1,subI}   = ERPallcond;
    allsubs_ERPtwocond_keeptrials{1,subI}   = ERPtwocond;
    allsubs_allfirstERPs_keeptrials{1,subI} = allfirstERPs;
    allsubs_alllastERPs_keeptrials{1,subI}  = alllastERPs;
   
    %% Plot the ERPs
    
    % ------------------------------------
    % a) plot over frontal sensors (diff vs easy) averaged ERPs
    cfg = [];
    cfg.channel = {'F1', 'F3', 'F5', 'F7', 'Fz', 'F2', 'F4', 'F6', 'F8'};
    figure; 
    ft_singleplotER(cfg, avgERPtwocond{1,1}, avgERPtwocond{1,2})
    title('forntal channels averaged ERPs')
    legend('easy', 'difficult')
    
    print(gcf, '-dpng', ['beads_analysis/figures/erps/beads_sub', num2str(subI), '_fig1_avrERP_frontal'])

    % b) plot over parietal sensors (diff vs easy) averaged ERPs
    cfg = [];
    cfg.channel = {'P1', 'P3', 'P5', 'P7', 'P9', 'Pz', 'P2', 'P4', 'P6', 'P8', 'P10'};
    figure; 
    ft_singleplotER(cfg, avgERPtwocond{1,1}, avgERPtwocond{1,2})
    title('parietal channels averaged ERPs')
    legend('easy', 'difficult')
    
    print(gcf, '-dpng', ['beads_analysis/figures/erps/beads_sub', num2str(subI), '_fig2_avrERP_parietal'])

    % -------------------------------------------
    % c) plot [all draws - last] vs last draw for easy cond (over frontal
    % electrodes)
    cfg = [];
    cfg.channel = {'F1', 'F3', 'F5', 'F7', 'Fz', 'F2', 'F4', 'F6', 'F8'};
    figure; 
    ft_singleplotER(cfg, allButLastERPs{1,1}, lastERPs{1,1})
    title('all draws but last vs last draw easy-frontal')
    legend('allButLast', 'last')
    
    % save figures
    print(gcf, '-dpng', ['beads_analysis/figures/erps/beads_sub', num2str(subI), '_frontal_easy_allButLastERPs'])

    % d) plot [all draws - last] vs last draw for easy cond (over parietal
    % electrodes)
    cfg = [];
    cfg.channel = {'P1', 'P3', 'P5', 'P7', 'P9', 'Pz', 'P2', 'P4', 'P6', 'P8', 'P10'};
    figure; 
    ft_singleplotER(cfg, allButLastERPs{1,1}, lastERPs{1,1})
    title('all draws but last vs last draw easy-parietal')
    legend('allButLast', 'last')
    
    print(gcf, '-dpng', ['beads_analysis/figures/erps/beads_sub', num2str(subI), '_parietal_easy_allButLastERPs'])
    
    % e) plot [all draws - last] vs last draw for diff cond (over frontal
    % electrodes)
    cfg = [];
    cfg.channel = {'F1', 'F3', 'F5', 'F7', 'Fz', 'F2', 'F4', 'F6', 'F8'};
    figure; 
    ft_singleplotER(cfg, allButLastERPs{1,2}, lastERPs{1,2})
    title('all draws but last vs last draw diff-frontal')
    legend('allButLast', 'last')
    
    print(gcf, '-dpng', ['beads_analysis/figures/erps/beads_sub', num2str(subI), '_frontal_diff_allButLastERPs'])
    
    % f) plot [all draws - last] vs last draw for diff cond (over parietal
    % electrodes)
    cfg = [];
    cfg.channel = {'P1', 'P3', 'P5', 'P7', 'P9', 'Pz', 'P2', 'P4', 'P6', 'P8', 'P10'};
    figure; 
    ft_singleplotER(cfg, allButLastERPs{1,2}, lastERPs{1,2})
    title('all draws but last vs last draw diff-parietal')
    legend('allButLast', 'last')
    
    print(gcf, '-dpng', ['beads_analysis/figures/erps/beads_sub', num2str(subI), '_parietal_diff_allButLastERPs'])
    
    % ----------------------------------------------------------
    % plot ERPs on interactive mode (plot all but last and last for easy cond)
    cfg = [];
    cfg.layout  = 'biosemi64.lay';
    cfg.interactive = 'yes';
    figure; ft_multiplotER(cfg, allButLastERPs{1,1}, lastERPs{1,1})
    
    % plot ERPs on interactive mode (plot all but last and last draws for diff cond)
    cfg = [];
    cfg.layout  = 'biosemi64.lay';
    cfg.interactive = 'yes';
    figure; ft_multiplotER(cfg, allButLastERPs{1,2}, lastERPs{1,2})
    
    % -----------------------------------------------------------------
    % maybe plot main ERPs in topoplots?
    % ALL but last, EASY
    cfg = [];
    cfg.xlim = [0.3 0.5];
    % cfg.zlim = [0 6e-14]; % this command messes up with the plot for some reason 
    cfg.layout  = 'biosemi64.lay';
    cfg.parameter = 'avg';
    cfg.interactive = 'yes';
    figure; ft_topoplotER(cfg, allButLastERPs{1,1}); colorbar
    
    % save
    print(gcf, '-dpng', ['beads_analysis/figures/erps/beads_sub', num2str(subI), '_allbutlast_easy_topo'])
    
    % Last, Easy
    cfg = [];
    cfg.xlim = [0.3 0.5];
    cfg.layout  = 'biosemi64.lay';
    cfg.parameter = 'avg';
    cfg.interactive = 'yes';
    figure; ft_topoplotER(cfg, lastERPs{1,1}); colorbar
    
    print(gcf, '-dpng', ['beads_analysis/figures/erps/beads_sub', num2str(subI), '_last_easy_topo'])
    
    % ALL but last, DIFF
    cfg = [];
    cfg.xlim = [0.3 0.5];
    cfg.layout  = 'biosemi64.lay';
    cfg.parameter = 'avg';
    cfg.interactive = 'yes';
    figure; ft_topoplotER(cfg, allButLastERPs{1,2}); colorbar
    
    % save
    print(gcf, '-dpng', ['beads_analysis/figures/erps/beads_sub', num2str(subI), '_allbutlast_diff_topo'])
    
    % Last, DIFF
    cfg = [];
    cfg.xlim = [0.3 0.5];
    cfg.layout  = 'biosemi64.lay';
    cfg.parameter = 'avg';
    cfg.interactive = 'yes';
    figure; ft_topoplotER(cfg, lastERPs{1,2}); colorbar
    
    % save
    print(gcf, '-dpng', ['beads_analysis/figures/erps/beads_sub', num2str(subI), '_last_diff_topo'])
    
  
    %% Plot the differences (contrasts) using avg structures for each subject 

    % partietal sensors between easy and diff trials over frontal &
    % parietal sites 
    cfg                     = [];
    cfg.operation           = 'x2-x1';
    cfg.parameter           = 'avg';
    difference_condwave     = ft_math(cfg, avgERPtwocond{1,1}, avgERPtwocond{1,2});
    
    % plot over frontal sensors
    cfg = [];
    cfg.channel = {'F1', 'F3', 'F5', 'F7', 'Fz', 'F2', 'F4', 'F6', 'F8'};
    figure; 
    ft_singleplotER(cfg, difference_condwave)
    title('forntal channels difference in averaged ERPs')
    
    print(gcf, '-dpng', ['beads_analysis/figures/averages/beads_sub', num2str(subI), '_fig1_differenceERPavg_twoconds'])
    
    % plot over parietal sensors 
    cfg = [];
    cfg.channel = {'P1', 'P3', 'P5', 'P7', 'P9', 'Pz', 'P2', 'P4', 'P6', 'P8', 'P10'};
    figure; 
    ft_singleplotER(cfg, difference_condwave)
    title('parietal channels difference in averaged ERPs')
    
    print(gcf, '-dpng', ['beads_analysis/figures/averages/beads_sub', num2str(subI), '_fig2_differenceERPavg_twoconds'])
    
    % show difference between conditions in time as a movie
    figure
    cfg        = [];
    cfg.layout = 'biosemi64.lay';
    ft_movieplotER(cfg, difference_condwave); colorbar

    % -------------------------------------------------------
    % Plot the difference between all-but-last and last draws (for easy and diff trials) over frontal &
    % partietal sensors
    
    % first calculate the differences for all-but-last vs last draws ERPs
    % for easy and diff trials
    cfg                     = [];
    cfg.operation           = 'x2-x1';
    cfg.parameter           = 'avg';
    difference_wave_one     = ft_math(cfg, allButLastERPs{1,1}, lastERPs{1,1}); % for easy trials 
    difference_wave_two     = ft_math(cfg, allButLastERPs{1,2}, lastERPs{1,2}); % for diff trials 
    
    % plot over frontal sensors
    cfg = [];
    cfg.channel = {'F1', 'F3', 'F5', 'F7', 'Fz', 'F2', 'F4', 'F6', 'F8'};
    figure; 
    ft_singleplotER(cfg, difference_wave_one, difference_wave_two)
    title('forntal channels difference in all-but-last vs last ERPs')
    legend('easy', 'difficult')
    
    print(gcf, '-dpng', ['beads_analysis/figures/averages/beads_sub', num2str(subI), '_fig1_differenceERPallbutlast_vs_last'])
    
    % plot over parietal sensors 
    cfg = [];
    cfg.channel = {'P1', 'P3', 'P5', 'P7', 'P9', 'Pz', 'P2', 'P4', 'P6', 'P8', 'P10'};
    figure; 
    ft_singleplotER(cfg, difference_wave_one, difference_wave_two)
    title('parietal channels difference in all-but-last vs last ERPs')
    legend('easy', 'difficult')
    
    print(gcf, '-dpng', ['beads_analysis/figures/averages/beads_sub', num2str(subI), '_fig2_differenceERPallbutlast_vs_last'])
    
    % plot diferences as a movie
    % differences between all-but-last vs last draw ERPs for easy cond
    figure
    cfg        = [];
    cfg.layout = 'biosemi64.lay';
    ft_movieplotER(cfg, difference_wave_one); colorbar
    
    % differences between all-but-last vs last draw ERPs for diff cond
    figure
    cfg        = [];
    cfg.layout = 'biosemi64.lay';
    ft_movieplotER(cfg, difference_wave_two); colorbar
    
    
    
    %% Time-Frequency Representation Analysis (TFR)
    
    % For TFR analysis I will need to use the "partdata" structure (the
    % re-referenced data with the non-eeg channels removed) 
    
    % I will use Morlet Wavelets 
    % One crucial parameter is width, as it determines the width of the
    % wavelets in number of cycles. With width set to 7: 
    % spectral bandwidth of frequency f:f/width*2, so at 30Hz freq: 30/7*2 = 8.6Hz 
    % wavelet duration: width/f/pi, in this case 7/30/pi=0.074s (or 74ms)
    
    % split data into 2 groups of different categories:
    % Group 1: condition group: easy vs diff
    % Group 2: draws group: all But Last vs Last draws
    
    % Group 1: run TFR for easy vs diff 
    for blocki = 1:blocks
        
        % load structure 
        disp(['loading beads_analysis/prepro/beads_preproc_sub_', num2str(subI), '_block_', num2str(blocki)])
        load(['beads_analysis/prepro/beads_preproc_sub_', num2str(subI), '_block_', num2str(blocki)], 'data')

        easy_data(blocki)   = data.easydata;
        diff_data(blocki)   = data.diffdata;

    end
    
    cfg             = [];
    all_easydata    = ft_appenddata(cfg, easy_data(1), easy_data(2), easy_data(3), easy_data(4));
    all_diffdata    = ft_appenddata(cfg, diff_data(1), diff_data(2), diff_data(3), diff_data(4));
    
    %% Run TFR and keep trials (will be used for between-trials stats)
    
    % ---------------------------
    % run TFR for Group 1 conditions
    cfg = [];
    cfg.channel    = 'all';
    cfg.method     = 'wavelet';
    cfg.width      = 7; % width definition is important 
    cfg.pad        = 10;
    cfg.output     = 'pow';
    cfg.foi        = 1:1:30;
    cfg.toi        = 'all'; % for computation efficiency use 'all'
    cfg.keeptrials = 'yes';
    TFRwave_easy   = ft_freqanalysis(cfg, all_easydata);
    TFRwave_diff   = ft_freqanalysis(cfg, all_diffdata);
    
    % save mat files in all subs cell 
    allsub_TFReasy{1,subI} = TFRwave_easy;
    allsub_TFRdiff{1,subI} = TFRwave_diff;
     
    % ---------------------------
    % run TFR for Group 2 conditions (all-but-last draws vs last darw)
    cfg = [];
    cfg.channel    = 'all';
    cfg.method     = 'wavelet';
    cfg.width      = 7; % width definition is important 
    cfg.pad        = 10;
    cfg.output     = 'pow';
    cfg.foi        = 1:1:30;
    cfg.toi        = 'all'; % for computation efficiency use 'all'
    cfg.keeptrials = 'yes';
    TFRwave_first  = ft_freqanalysis(cfg, all_firstdata);
    TFRwave_last   = ft_freqanalysis(cfg, all_lastdata);
    
    % save mat files in all subs cell 
    allsub_TFRfirst{1,subI}     = TFRwave_first;
    allsub_TFRlast{1,subI}      = TFRwave_last;
    
    % ---------------------------
    % run TFR for Group 2 conditions (easy all-but-last and last darws) and
    % (diff all-but-last and last darws)
    for icond = 1:totalconds
        
        cfg                         = [];
        cfg.channel                 = 'all';
        cfg.method                  = 'wavelet';
        cfg.width                   = 7; % width definition is important 
        cfg.pad                     = 10;
        cfg.output                  = 'pow';
        cfg.foi                     = 1:1:30;
        cfg.toi                     = 'all'; % for computation efficiency use 'all'
        cfg.trials                  = find(all_firstdata.trialinfo(:,4) == icond);
        cfg.keeptrials              = 'yes';
        
        TFRwave_first_2cond{icond}  = ft_freqanalysis(cfg, all_firstdata);
        
        cfg = [];
        cfg.channel                 = 'all';
        cfg.method                  = 'wavelet';
        cfg.width                   = 7; % width definition is important 
        cfg.pad                     = 10;
        cfg.output                  = 'pow';
        cfg.foi                     = 1:1:30;
        cfg.toi                     = 'all'; % for computation efficiency use 'all'
        cfg.trials                  = find(all_lastdata.trialinfo(:,4) == icond);
        cfg.keeptrials              = 'yes';
        
        TFRwave_last_2cond{icond}   = ft_freqanalysis(cfg, all_lastdata);

    end
    
    % save files in all subs cell 
    allsub_TFRfirst_2cond{1,subI}   = TFRwave_first_2cond;
    allsub_TFRlast_2cond{1,subI}    = TFRwave_last_2cond;
    
    
    %% Run TFR averaged for Within-Subjects analysis 
    
    % this analysis will be the same as the first, only we will not keep
    % trials 
    
    % ---------------------------
    % run TFR for Group 1 conditions
    cfg                 = [];
    cfg.channel         = 'all';
    cfg.method          = 'wavelet';
    cfg.width           = 7; % width definition is important 
    cfg.pad             = 10;
    cfg.output          = 'pow';
    cfg.foi             = 1:1:30;
    cfg.toi             = 'all'; % for computation efficiency use 'all'
    TFRwave_easy_avg    = ft_freqanalysis(cfg, all_easydata);
    TFRwave_diff_avg    = ft_freqanalysis(cfg, all_diffdata);
    
    % save mat files in all subs cell 
    allsub_TFReasy_avg{1,subI} = TFRwave_easy_avg;
    allsub_TFRdiff_avg{1,subI} = TFRwave_diff_avg;

    cfg = [];
    cfg.baseline     = [-0.2 0];
    cfg.baselinetype = 'absolute';
    cfg.marker       = 'on';
    cfg.showlabels   = 'yes';
    cfg.layout       = 'biosemi64.lay';
    figure; ft_multiplotTFR(cfg, TFRwave_easy_avg); ft_multiplotTFR(cfg, TFRwave_diff_avg); colorbar;
    
    % single plot 
    cfg = [];
    cfg.baseline     = [-0.2 0];
    cfg.baselinetype = 'absolute';
    cfg.marker       = 'on';
    cfg.maskstyle    = 'saturation';
    cfg.channel      = {'P1', 'P3', 'P5', 'P7', 'P9', 'Pz', 'P2', 'P4', 'P6', 'P8', 'P10'};
    cfg.layout       = 'biosemi64.lay';
    
    figure; ft_singleplotTFR(cfg, TFRwave_easy_avg); colorbar;
    print(gcf, '-dpng', ['beads_analysis/figures/tfr/beads_sub', num2str(subI), '_tfr_easy'])
    
    figure; ft_singleplotTFR(cfg, TFRwave_diff_avg); colorbar;
    print(gcf, '-dpng', ['beads_analysis/figures/tfr/beads_sub', num2str(subI), '_tfr_diff'])
    
    % ---------------------------
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
    cfg.maskstyle        = 'saturation';
    cfg.channel      = {'P1', 'P3', 'P5', 'P7', 'P9', 'Pz', 'P2', 'P4', 'P6', 'P8', 'P10'};
    cfg.layout       = 'biosemi64.lay';
    
    figure; ft_singleplotTFR(cfg, TFRwave_first_avg); colorbar;
    print(gcf, '-dpng', ['beads_analysis/figures/tfr/beads_sub', num2str(subI), '_tfr_all_but_last'])
    
    figure; ft_singleplotTFR(cfg, TFRwave_last_avg); colorbar;
    print(gcf, '-dpng', ['beads_analysis/figures/tfr/beads_sub', num2str(subI), '_tfr_last'])
    
    % ---------------------------
    % run TFR for Group 2 conditions (easy all-but-last and last darws) and
    % (diff all-but-last and last darws)
    for icond = 1:totalconds
        
        cfg                             = [];
        cfg.channel                     = 'all';
        cfg.method                      = 'wavelet';
        cfg.width                       = 7; % width definition is important 
        cfg.pad                         = 10;
        cfg.output                      = 'pow';
        cfg.foi                         = 1:1:30;
        cfg.toi                         = 'all'; % for computation efficiency use 'all'
        cfg.trials                      = find(all_firstdata.trialinfo(:,4) == icond);
        
        TFRwave_first_2cond_avg{icond}  = ft_freqanalysis(cfg, all_firstdata);
        
        cfg                             = [];
        cfg.channel                     = 'all';
        cfg.method                      = 'wavelet';
        cfg.width                       = 7; % width definition is important 
        cfg.pad                         = 10;
        cfg.output                      = 'pow';
        cfg.foi                         = 1:1:30;
        cfg.toi                         = 'all'; % for computation efficiency use 'all'
        cfg.trials                      = find(all_lastdata.trialinfo(:,4) == icond);
        
        TFRwave_last_2cond_avg{icond}   = ft_freqanalysis(cfg, all_lastdata);

    end
    
    % save files in all subs cell 
    allsub_TFRfirst_2cond_avg{1,subI}   = TFRwave_first_2cond_avg;
    allsub_TFRlast_2cond_avg{1,subI}    = TFRwave_last_2cond_avg;
   
end % end of subject loop

%% Save all sub cell files for later subject-level & group-level analyses

% save all subjects ERP analyes 
save('beads_analysis/erps/beads_allsubsAllCondsERPs', 'allsubERPallcond')
save('beads_analysis/erps/beads_allsubsTwoCondsERPs', 'allsubERPtwocond')
save('beads_analysis/erps/beads_allsubsAllButLastERPs', 'allsubfirstERPs')
save('beads_analysis/erps/beads_allsubsLastERPs', 'allsublastERPs')

save('beads_analysis/erps/beads_allsubs_keeptrialsAllERPs', 'allsubs_ERPallcond_keeptrials')
save('beads_analysis/erps/beads_allsubs_keeptrialsTwoERPs', 'allsubs_ERPtwocond_keeptrials')
save('beads_analysis/erps/beads_allsubs_keeptrialsFirstERPs', 'allsubs_allfirstERPs_keeptrials')
save('beads_analysis/erps/beads_allsubs_keeptrialsLastERPs', 'allsubs_alllastERPs_keeptrials')

% save all subjects TFR analyes - with the full trials
save('beads_analysis/tfr/beads_allsubs_tfrEasy', 'allsub_TFReasy', '-v7.3')
save('beads_analysis/tfr/beads_allsubs_tfrDiff', 'allsub_TFRdiff', '-v7.3')
save('beads_analysis/tfr/beads_allsubs_tfrFirst', 'allsub_TFRfirst', '-v7.3')
save('beads_analysis/tfr/beads_allsubs_tfrLast', 'allsub_TFRlast', '-v7.3')
save('beads_analysis/tfr/beads_allsubs_tfrFirst_2cond', 'allsub_TFRfirst_2cond', '-v7.3')
save('beads_analysis/tfr/beads_allsubs_tfrLast_2cond', 'allsub_TFRlast_2cond', '-v7.3')

% save all subject TFR analyses - averages
save('beads_analysis/tfr/beads_allsubs_tfrEasy_avg', 'allsub_TFReasy_avg', '-v7.3')
save('beads_analysis/tfr/beads_allsubs_tfrDiff_avg', 'allsub_TFRdiff_avg', '-v7.3')
save('beads_analysis/tfr/beads_allsubs_tfrFirst_avg', 'allsub_TFRfirst_avg', '-v7.3')
save('beads_analysis/tfr/beads_allsubs_tfrLast_avg', 'allsub_TFRlast_avg', '-v7.3')
save('beads_analysis/tfr/beads_allsubs_tfrFirst_2cond_avg', 'allsub_TFRfirst_2cond_avg', '-v7.3')
save('beads_analysis/tfr/beads_allsubs_tfrLast_2cond_avg', 'allsub_TFRlast_2cond_avg', '-v7.3')

%% Compute TFR grand averages for all datasets/conditions

% calculate grand average for each condition with the individual data set
% to yes
cfg                 = [];
cfg.channel         = 'all';
cfg.toilim          = 'all';
cfg.foilim          = 'all';
cfg.keepindividual  = 'yes';
cfg.parameter       = 'powspctrm';

ga_TFReasy_ind      = ft_freqgrandaverage(cfg, allsub_TFReasy_avg{:});
ga_TFRdiff_ind      = ft_freqgrandaverage(cfg, allsub_TFRdiff_avg{:});
ga_TFRfirst_ind     = ft_freqgrandaverage(cfg, allsub_TFRfirst_avg{:});
ga_TFRlast_ind      = ft_freqgrandaverage(cfg, allsub_TFRlast_avg{:});

% split easy and difficult structures for each subject
for i = 1:length(allsub_TFRlast_avg)
        
    tfr_tmp_first_easy{1,i} = allsub_TFRfirst_2cond_avg{1,i}{1,1};
    tfr_tmp_last_easy{1,i}  = allsub_TFRlast_2cond_avg{1,i}{1,1};
    
    tfr_tmp_first_diff{1,i} = allsub_TFRfirst_2cond_avg{1,i}{1,2};
    tfr_tmp_last_diff{1,i}  = allsub_TFRlast_2cond_avg{1,i}{1,2};

end

ga_TFRfirst_easy_ind        = ft_freqgrandaverage(cfg, tfr_tmp_first_easy{:});
ga_TFRlast_easy_ind         = ft_freqgrandaverage(cfg, tfr_tmp_last_easy{:});
ga_TFRfirst_diff_ind        = ft_freqgrandaverage(cfg, tfr_tmp_first_diff{:});
ga_TFRlast_diff_ind         = ft_freqgrandaverage(cfg, tfr_tmp_last_diff{:});

clear tfr_tmp_first_easy tfr_tmp_last_easy tfr_tmp_first_diff tfr_tmp_last_diff

% calculate grand average for each condition 
cfg                 = [];
cfg.channel         = 'all';
cfg.toilim          = 'all';
cfg.foilim          = 'all';
cfg.parameter       = 'powspctrm';

ga_TFReasy      = ft_freqgrandaverage(cfg, allsub_TFReasy_avg{:});
ga_TFRdiff      = ft_freqgrandaverage(cfg, allsub_TFRdiff_avg{:});
ga_TFRfirst     = ft_freqgrandaverage(cfg, allsub_TFRfirst_avg{:});
ga_TFRlast     = ft_freqgrandaverage(cfg, allsub_TFRlast_avg{:});

% split easy and difficult structures for each subject
for i = 1:length(allsub_TFRlast)
        
    tfr_tmp_first_easy{1,i} = allsub_TFRfirst_2cond_avg{1,i}{1,1};
    tfr_tmp_last_easy{1,i}  = allsub_TFRlast_2cond_avg{1,i}{1,1};
    
    tfr_tmp_first_diff{1,i} = allsub_TFRfirst_2cond_avg{1,i}{1,2};
    tfr_tmp_last_diff{1,i}  = allsub_TFRlast_2cond_avg{1,i}{1,2};

end

ga_TFRfirst_easy        = ft_freqgrandaverage(cfg, tfr_tmp_first_easy{:});
ga_TFRlast_easy         = ft_freqgrandaverage(cfg, tfr_tmp_last_easy{:});
ga_TFRfirst_diff        = ft_freqgrandaverage(cfg, tfr_tmp_first_diff{:});
ga_TFRlast_diff         = ft_freqgrandaverage(cfg, tfr_tmp_last_diff{:});

% save TFR grand averages (with individual trials) 
save('beads_analysis/tfr/beads_ga_TFReasy_ind', 'ga_TFReasy_ind', '-v7.3')
save('beads_analysis/tfr/beads_ga_TFRdiff_ind', 'ga_TFRdiff_ind', '-v7.3')
save('beads_analysis/tfr/beads_ga_TFRfirst_ind', 'ga_TFRfirst_ind', '-v7.3')
save('beads_analysis/tfr/beads_ga_TFRlast_ind', 'ga_TFRlast_ind', '-v7.3')
save('beads_analysis/tfr/beads_ga_TFRfirst_easy_ind', 'ga_TFRfirst_easy_ind', '-v7.3')
save('beads_analysis/tfr/beads_ga_TFRlast_easy_ind', 'ga_TFRlast_easy_ind', '-v7.3')
save('beads_analysis/tfr/beads_ga_TFRfirst_diff_ind', 'ga_TFRfirst_diff_ind', '-v7.3')
save('beads_analysis/tfr/beads_ga_TFRlast_diff_ind', 'ga_TFRlast_diff_ind', '-v7.3')

% save TFR grand averages (without individual trials) 
save('beads_analysis/tfr/beads_ga_TFReasy', 'ga_TFReasy', '-v7.3')
save('beads_analysis/tfr/beads_ga_TFRdiff', 'ga_TFRdiff', '-v7.3')
save('beads_analysis/tfr/beads_ga_TFRfirst', 'ga_TFRfirst', '-v7.3')
save('beads_analysis/tfr/beads_ga_TFRlast', 'ga_TFRlast', '-v7.3')
save('beads_analysis/tfr/beads_ga_TFRfirst_easy', 'ga_TFRfirst_easy', '-v7.3')
save('beads_analysis/tfr/beads_ga_TFRlast_easy', 'ga_TFRlast_easy', '-v7.3')
save('beads_analysis/tfr/beads_ga_TFRfirst_diff', 'ga_TFRfirst_diff', '-v7.3')
save('beads_analysis/tfr/beads_ga_TFRlast_diff', 'ga_TFRlast_diff', '-v7.3')


%% Compute Timelock Grand Averages for each dataset/condition 

% first split the ERP to seperate cells for easy vs difficult 
for i = 1:length(allsubfirstERPs)

    tmp_first_easy{1,i}     = allsubfirstERPs{1,i}{1,1};
    tmp_last_easy{1,i}      = allsublastERPs{1,i}{1,1};

    tmp_first_diff{1,i}     = allsubfirstERPs{1,i}{1,2};
    tmp_last_diff{1,i}      = allsublastERPs{1,i}{1,2};

    tmp_easy{1,i}           = allsubERPtwocond{1,i}{1,1};
    tmp_diff{1,i}           = allsubERPtwocond{1,i}{1,2};

end

% now calculate timelocked grand averages 
cfg             = [];
cfg.channel     = 'all';
cfg.latency     = 'all';
cfg.parameter   = 'avg';

gaERP_easy      = ft_timelockgrandaverage(cfg, tmp_easy{:});
gaERP_diff      = ft_timelockgrandaverage(cfg, tmp_diff{:});
gaERP_firsteasy = ft_timelockgrandaverage(cfg, tmp_first_easy{:});
gaERP_lasteasy  = ft_timelockgrandaverage(cfg, tmp_last_easy{:});
gaERP_firstdiff = ft_timelockgrandaverage(cfg, tmp_first_diff{:});
gaERP_lastdiff  = ft_timelockgrandaverage(cfg, tmp_last_diff{:});

% save timelock averages 
save('beads_analysis/erps/beads_gaERPeasy', 'gaERP_easy')
save('beads_analysis/erps/beads_gaERPdiff', 'gaERP_diff')
save('beads_analysis/erps/beads_gaERP_firsteasy', 'gaERP_firsteasy')
save('beads_analysis/erps/beads_gaERP_lasteasy', 'gaERP_lasteasy')
save('beads_analysis/erps/beads_gaERP_firstdiff', 'gaERP_firstdiff')
save('beads_analysis/erps/beads_gaERP_lastdiff', 'gaERP_lastdiff')

% load the ERP grand averages if needed 
% load('beads_analysis/erps/beads_gaERPeasy', 'gaERP_easy')
% load('beads_analysis/erps/beads_gaERPdiff', 'gaERP_diff')
% load('beads_analysis/erps/beads_gaERP_firsteasy', 'gaERP_firsteasy')
% load('beads_analysis/erps/beads_gaERP_lasteasy', 'gaERP_lasteasy')
% load('beads_analysis/erps/beads_gaERP_firstdiff', 'gaERP_firstdiff')
% load('beads_analysis/erps/beads_gaERP_lastdiff', 'gaERP_lastdiff')
    
%% Visualise grand average ERPs and compute GA contrasts

% plot grand averages for difficults vs easy ERPs over frontal sites 
cfg = [];
cfg.channel = {'F1', 'F3', 'F5', 'F7', 'Fz', 'F2', 'F4', 'F6', 'F8'};
figure; 
ft_singleplotER(cfg, gaERP_easy, gaERP_diff)
title('forntal channels GA ERPs')
legend('easy', 'difficult')

print(gcf, '-dpng', 'beads_analysis/figures/grand_averageERPs_plots/beads_GA_easy_vs_diff_ERP_frontal')

% ------------------------------

% plot grand averages for difficult vs easy ERPs over parietal sites 
cfg = [];
cfg.channel = {'P1', 'P3', 'P5', 'P7', 'P9', 'Pz', 'P2', 'P4', 'P6', 'P8', 'P10'};
figure; 
ft_singleplotER(cfg, gaERP_easy, gaERP_diff)
title('Parietal channels GA ERPs')
legend('easy', 'difficult')

print(gcf, '-dpng', 'beads_analysis/figures/grand_averageERPs_plots/beads_GA_easy_vs_diff_ERP_parietal')

% ---------------------------------
% plot grand averages for all-but-last vs last easy draws ERPs over frontal sites 
cfg = [];
cfg.channel = {'F1', 'F3', 'F5', 'F7', 'Fz', 'F2', 'F4', 'F6', 'F8'};
figure; 
ft_singleplotER(cfg, gaERP_firsteasy, gaERP_lasteasy)
title('GA ERPs all-draws-but-last vs last for easy-frontal')
legend('all-but-last', 'last')

print(gcf, '-dpng', 'beads_analysis/figures/grand_averageERPs_plots/beads_GA_first_vs_last_easy_ERP_frontal')

% ---------------------------------
% plot grand averages for all-but-last vs last easy draws ERPs over parietal sites 
cfg = [];
cfg.channel = {'P1', 'P3', 'P5', 'P7', 'P9', 'Pz', 'P2', 'P4', 'P6', 'P8', 'P10'};
figure; 
ft_singleplotER(cfg, gaERP_firsteasy, gaERP_lasteasy)
title('GA ERPs all-draws-but-last vs last for easy-parietal')
legend('all-but-last', 'last')

print(gcf, '-dpng', 'beads_analysis/figures/grand_averageERPs_plots/beads_GA_first_vs_last_easy_ERP_parietal')

% ---------------------------------
% plot grand averages for all-but-last vs last diff draws ERPs over frontal sites 
cfg = [];
cfg.channel = {'F1', 'F3', 'F5', 'F7', 'Fz', 'F2', 'F4', 'F6', 'F8'};
figure; 
ft_singleplotER(cfg, gaERP_firstdiff, gaERP_lastdiff)
title('GA ERPs all-draws-but-last vs last for difficult-frontal')
legend('all-but-last', 'last')

print(gcf, '-dpng', 'beads_analysis/figures/grand_averageERPs_plots/beads_GA_first_vs_last_difficult_ERP_frontal')

% ---------------------------------
% plot grand averages for all-but-last vs last easy draws ERPs over parietal sites 
cfg = [];
cfg.channel = {'P1', 'P3', 'P5', 'P7', 'P9', 'Pz', 'P2', 'P4', 'P6', 'P8', 'P10'};
figure; 
ft_singleplotER(cfg, gaERP_firstdiff, gaERP_lastdiff)
title('GA ERPs all-draws-but-last vs last for difficult-parietal')
legend('all-but-last', 'last')

print(gcf, '-dpng', 'beads_analysis/figures/grand_averageERPs_plots/beads_GA_first_vs_last_diff_ERP_parietal')

%% Compute GA differences for each codnition/dataset

% first compute the averages 
cfg                         = [];
cfg.operation               = 'x2-x1';
cfg.parameter               = 'avg';
difference_diffvseasy       = ft_math(cfg, gaERP_diff, gaERP_easy);
difference_lastvsfirst_easy = ft_math(cfg, gaERP_lasteasy, gaERP_firsteasy);
difference_lastvsfirst_diff = ft_math(cfg, gaERP_lastdiff, gaERP_firstdiff);

% compute differences for diff [last vs all-but-last] vs easy [last vs all-but-last]
cfg                         = [];
cfg.operation               = '(x2-x1) - (x4-x3)';
cfg.parameter               = 'avg';
difference_lastvsfirst      = ft_math(cfg, gaERP_lastdiff, gaERP_firstdiff, gaERP_lasteasy, gaERP_firsteasy);
    
% plot difference between easy vs difficult over frontal sites 
cfg = [];
cfg.channel = {'F1', 'F3', 'F5', 'F7', 'Fz', 'F2', 'F4', 'F6', 'F8'};
figure; 
ft_singleplotER(cfg, difference_diffvseasy)
title('GA difference easy vs difficult ERPs over frontal sites')

print(gcf, '-dpng', 'beads_analysis/figures/differences/beads_figGA_difference_easyvsdiff_frontal')

% plot difference between easy vs difficult over frontal sites 
cfg = [];
cfg.channel = {'P1', 'P3', 'P5', 'P7', 'P9', 'Pz', 'P2', 'P4', 'P6', 'P8', 'P10'};
figure; 
ft_singleplotER(cfg, difference_diffvseasy)
title('GA difference easy vs difficult ERPs over parietal sites')

print(gcf, '-dpng', 'beads_analysis/figures/differences/beads_figGA_difference_easyvsdiff_parietal')

% plot GA difference between last vs all-but-last ERP draws over frontal
% ellectrodes for easy vs difficult
cfg = [];
cfg.channel = {'F1', 'F3', 'F5', 'F7', 'Fz', 'F2', 'F4', 'F6', 'F8'};
figure; 
ft_singleplotER(cfg, difference_lastvsfirst_easy, difference_lastvsfirst_diff)
title('GA difference in all-but-last vs last ERPs for forntal sites')
legend('easy', 'difficult')

print(gcf, '-dpng', 'beads_analysis/figures/differences/beads_figGA_difference_lastvsfirst_frontal')
    
% plot GA difference between last vs all-but-last ERP draws over parietal
% ellectrodes for easy vs difficult
cfg = [];
cfg.channel = {'P1', 'P3', 'P5', 'P7', 'P9', 'Pz', 'P2', 'P4', 'P6', 'P8', 'P10'};
figure; 
ft_singleplotER(cfg, difference_lastvsfirst_easy, difference_lastvsfirst_diff)
title('GA difference in all-but-last vs last ERPs for parietal sites')
legend('easy', 'difficult')

print(gcf, '-dpng', 'beads_analysis/figures/differences/beads_figGA_difference_lastvsfirst_parietal')

% plot difference between easy vs difficult over frontal sites 
cfg = [];
cfg.channel = {'F1', 'F3', 'F5', 'F7', 'Fz', 'F2', 'F4', 'F6', 'F8'};
figure; 
ft_singleplotER(cfg, difference_lastvsfirst)
title('GA total difference easy vs difficult ERPs over frontal sites')

print(gcf, '-dpng', 'beads_analysis/figures/differences/beads_figGA_total_difference_easyvsdiff_frontal')

% plot difference between easy vs difficult over frontal sites 
cfg = [];
cfg.channel = {'P1', 'P3', 'P5', 'P7', 'P9', 'Pz', 'P2', 'P4', 'P6', 'P8', 'P10'};
figure; 
ft_singleplotER(cfg, difference_diffvseasy)
title('GA total difference easy vs difficult ERPs over parietal sites')

print(gcf, '-dpng', 'beads_analysis/figures/differences/beads_figGA_total_difference_easyvsdiff_parietal')

% plot total GA difference between last vs all-but-last as a movie 
figure
cfg        = [];
cfg.layout = 'biosemi64.lay';
ft_movieplotER(cfg, difference_diffvseasy); colorbar

%% Visualise grand average TFRs for all conditions/datasets

% load TFR GA if needed
load('beads_analysis/tfr/beads_ga_TFReasy', 'ga_TFReasy')
load('beads_analysis/tfr/beads_ga_TFRdiff', 'ga_TFRdiff')
load('beads_analysis/tfr/beads_ga_TFRfirst', 'ga_TFRfirst')
load('beads_analysis/tfr/beads_ga_TFRlast', 'ga_TFRlast')
load('beads_analysis/tfr/beads_ga_TFRfirst_easy', 'ga_TFRfirst_easy')
load('beads_analysis/tfr/beads_ga_TFRlast_easy', 'ga_TFRlast_easy')
load('beads_analysis/tfr/beads_ga_TFRfirst_diff', 'ga_TFRfirst_diff')
load('beads_analysis/tfr/beads_ga_TFRlast_diff', 'ga_TFRlast_diff')

% plot GA TFR for easy and difficult epochs over all sensors 
cfg                 = [];
cfg.baseline        = [-0.2 0];
cfg.baselinetype    = 'absolute';
cfg.marker          = 'on';
cfg.showlabels      = 'yes';
cfg.layout          = 'biosemi64.lay';
figure; ft_multiplotTFR(cfg, ga_TFReasy); colorbar; ft_multiplotTFR(cfg, ga_TFRdiff); colorbar;

% single plot of TFR easy and difficult over parietal sites
cfg = [];
cfg.baseline     = [-0.2 0];
cfg.baselinetype = 'absolute';
cfg.marker       = 'on';
cfg.maskstyle    = 'saturation';
cfg.channel      = {'P1', 'P3', 'P5', 'P7', 'P9', 'Pz', 'P2', 'P4', 'P6', 'P8', 'P10'};
cfg.layout       = 'biosemi64.lay';

figure; ft_singleplotTFR(cfg, ga_TFReasy); colorbar;
print(gcf, '-dpng', 'beads_analysis/figures/grand_averageTFRs_plots/beads_tfr_GA_easy')

figure; ft_singleplotTFR(cfg, ga_TFRdiff); colorbar;
print(gcf, '-dpng', 'beads_analysis/figures/grand_averageTFRs_plots/beads_tfr_GA_diff')

% plot GA TFR for easy all-but-last and last epochs over all sensors 
cfg                 = [];
cfg.baseline        = [-0.2 0];
cfg.baselinetype    = 'absolute';
cfg.marker          = 'on';
cfg.showlabels      = 'yes';
cfg.layout          = 'biosemi64.lay';
figure; ft_multiplotTFR(cfg, ga_TFRfirst_easy); colorbar; ft_multiplotTFR(cfg, ga_TFRlast_easy); colorbar;

% single plot of TFR easy all-but-last and last over parietal sites
cfg = [];
cfg.baseline     = [-0.2 0];
cfg.baselinetype = 'absolute';
cfg.marker       = 'on';
cfg.maskstyle    = 'saturation';
cfg.channel      = {'P1', 'P3', 'P5', 'P7', 'P9', 'Pz', 'P2', 'P4', 'P6', 'P8', 'P10'};
cfg.layout       = 'biosemi64.lay';

figure; ft_singleplotTFR(cfg, ga_TFRfirst_easy); colorbar;
print(gcf, '-dpng', 'beads_analysis/figures/grand_averageTFRs_plots/beads_tfr_GA_first_easy')

figure; ft_singleplotTFR(cfg, ga_TFRlast_easy); colorbar;
print(gcf, '-dpng', 'beads_analysis/figures/grand_averageTFRs_plots/beads_tfr_GA_last_easy')

% plot GA TFR for difficult all-but-last and last epochs over all sensors 
cfg                 = [];
cfg.baseline        = [-0.2 0];
cfg.baselinetype    = 'absolute';
cfg.marker          = 'on';
cfg.showlabels      = 'yes';
cfg.layout          = 'biosemi64.lay';
figure; ft_multiplotTFR(cfg, ga_TFRfirst_diff); colorbar; ft_multiplotTFR(cfg, ga_TFRlast_diff); colorbar;

% single plot of TFR easy all-but-last and last over parietal sites
cfg = [];
cfg.baseline     = [-0.2 0];
cfg.baselinetype = 'absolute';
cfg.marker       = 'on';
cfg.maskstyle    = 'saturation';
cfg.channel      = {'P1', 'P3', 'P5', 'P7', 'P9', 'Pz', 'P2', 'P4', 'P6', 'P8', 'P10'};
cfg.layout       = 'biosemi64.lay';

figure; ft_singleplotTFR(cfg, ga_TFRfirst_diff); colorbar;
print(gcf, '-dpng', 'beads_analysis/figures/grand_averageTFRs_plots/beads_tfr_GA_first_diff')

figure; ft_singleplotTFR(cfg, ga_TFRlast_diff); colorbar;
print(gcf, '-dpng', 'beads_analysis/figures/grand_averageTFRs_plots/beads_tfr_GA_last_diff')

% plot GA TFR for total all-but-last and last epochs over all sensors 
cfg                 = [];
cfg.baseline        = [-0.2 0];
cfg.baselinetype    = 'absolute';
cfg.marker          = 'on';
cfg.showlabels      = 'yes';
cfg.layout          = 'biosemi64.lay';
figure; ft_multiplotTFR(cfg, ga_TFRfirst); colorbar; ft_multiplotTFR(cfg, ga_TFRlast); colorbar;

% single plot of TFR total all-but-last and last over parietal sites
cfg = [];
cfg.baseline     = [-0.2 0];
cfg.baselinetype = 'absolute';
cfg.marker       = 'on';
cfg.maskstyle    = 'saturation';
cfg.channel      = {'P1', 'P3', 'P5', 'P7', 'P9', 'Pz', 'P2', 'P4', 'P6', 'P8', 'P10'};
cfg.layout       = 'biosemi64.lay';

figure; ft_singleplotTFR(cfg, ga_TFRfirst); colorbar;
print(gcf, '-dpng', 'beads_analysis/figures/grand_averageTFRs_plots/beads_tfr_GA_first')

figure; ft_singleplotTFR(cfg, ga_TFRlast); colorbar;
print(gcf, '-dpng', 'beads_analysis/figures/grand_averageTFRs_plots/beads_tfr_GA_last')


%% Run ERP statistics (between-trials/first level)

% 
load('beads_analysis/erps/beads_allsubs_keeptrialsAllERPs', 'allsubs_ERPallcond_keeptrials')
load('beads_analysis/erps/beads_allsubs_keeptrialsTwoERPs', 'allsubs_ERPtwocond_keeptrials')
load('beads_analysis/erps/beads_allsubs_keeptrialsFirstERPs', 'allsubs_allfirstERPs_keeptrials')
load('beads_analysis/erps/beads_allsubs_keeptrialsLastERPs', 'allsubs_alllastERPs_keeptrials')

% ------------------------------
% Run statitics at the trial-level for each subject first 
% Claster-based permutation tests will be used to look at differences
% between conditions: (easy vs diff) and draws: (all-but-last vs last darw)

% extract current subject mat file from a given cell loaded above and run
% analysis
for sub = 1:nsubs 
    
    % load the original easy and difficult data for current sub 
    for blocki = 1:blocks
        
        % load structure 
        disp(['loading beads_analysis/prepro/beads_preproc_sub_', num2str(sub), '_block_', num2str(blocki)])
        load(['beads_analysis/prepro/beads_preproc_sub_', num2str(sub), '_block_', num2str(blocki)], 'data')

        easy_data(blocki)   = data.easydata;
        diff_data(blocki)   = data.diffdata;

    end
    
    cfg                     = [];
    all_easydata            = ft_appenddata(cfg, easy_data(1), easy_data(2), easy_data(3), easy_data(4));
    all_diffdata            = ft_appenddata(cfg, diff_data(1), diff_data(2), diff_data(3), diff_data(4));
    
    % run cluster-based permutations using the monte-carlo method for easy
    % vs difficult trials/epochs for each subject seperately 
    thissub_twocondERPs     = allsubs_ERPtwocond_keeptrials{1,sub};
    
    condeasy                = thissub_twocondERPs{1,1};
    conddiff                = thissub_twocondERPs{1,2};
    
    % create the neighbours structure which will be critical to compute
    % stats
    cfg_neighb              = [];
    cfg_neighb.layout       = 'biosemi64.lay'; %in meters
    cfg_neighb.method       = 'distance';
    neighbours              = ft_prepare_neighbours(cfg_neighb, all_easydata);
    
    cfg                     = [];
    cfg.neighbours          = neighbours;   % the neighbours specify for each sensor with
                                            % which other sensors it can form clusters
    cfg.channel             = 'all';        % cell-array with selected channel labels
    cfg.latency             = [0.3 0.6];    % time interval over which the experimental
                                            % conditions must be compared (in seconds)

    % permutation tests
    cfg.method              = 'montecarlo';     % use the Monte Carlo Method to calculate the significance probability
    cfg.statistic           = 'indepsamplesT';  % use the independent samples T-statistic as a measure to
                                                % evaluate the effect at the sample level
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

    n_easy                  = size(condeasy.trial,1);
    n_diff                  = size(conddiff.trial,1);
    
    cfg.design              = [ones(1,n_diff), ones(1,n_easy)*2]; % design matrix
    cfg.ivar                = 1; % number or list with indices indicating the independent variable(s)

    % run stats 
    [stat]               = ft_timelockstatistics(cfg, conddiff, condeasy);
    
    %%% -------------------------------------------------------------- %%%
    % Make a vector of all p-values associated with the clusters from ft_timelockstatistics.
    pos_cluster_pvals = [stat.posclusters(:).prob];

    % Then, find which clusters are deemed interesting to visualize, here we use a cutoff criterion based on the
    % cluster-associated p-value, and take a 5% two-sided cutoff (i.e. 0.025 for the positive and negative clusters,
    % respectively
    pos_clust = find(pos_cluster_pvals < 0.025);
    pos       = ismember(stat.posclusterslabelmat, pos_clust);


    % and now for the negative clusters...
    neg_cluster_pvals = [stat.negclusters(:).prob];
    neg_clust         = find(neg_cluster_pvals < 0.025);
    neg               = ismember(stat.negclusterslabelmat, neg_clust);

    % if any positive cluster survived the cutoff.. loop over all sig positive clusters
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
            print(gcf, '-dpng', ['beads_analysis/figures/stats_erps/beads_fig1_STAT_pos', num2str(i), '_subject', num2str(sub)])
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
            print(gcf, '-dpng', ['beads_analysis/figures/stats_erps/beads_fig1_STAT_neg', num2str(i), '_subject', num2str(sub)])
        end
    end
    %%% -------------------------------------------------------------- %%%
    
    % save stats 
    save(['beads_analysis/stats/beads_sub', num2str(sub), '_stats_1stlevel_diff_vs_easy'], 'stat')
    
    % try ploting on clusterplots for time-windows, more informative than the plots above!!
    % negative & positive clusters
    % find max and min t-values first 
    if ~isempty(pos_clust) | ~isempty(neg_clust) % if there are any positive or negative clusters 
        maxval                  = max(stat.stat, [], 'all');
        minval                  = min(stat.stat, [], 'all');

        cfg                     = [];
        cfg.zlim                = [minval maxval]; % T-values
        cfg.alpha               = 0.025;
        cfg.highlightcolorpos   = [0 0 0];
        cfg.highlightcolorneg   = [1 1 0];
        cfg.saveaspng           = ['beads_analysis/figures/stats_erps/sub', num2str(sub), '_1stlevel_stats_diff_vs_easy'];
        cfg.toi                 = 0.3:0.01:0.6; % times of interest
        cfg.layout              = 'biosemi64.lay';
        ft_clusterplot(cfg,stat); 
    end

    clear stat neighbours cfg_neighb cfg
    
    %% Compute statistics for all-but-last vs last draws for easy trials 
    
    % load the original data for this subject
    for block = 1:blocks
        
        % if data structure is already in workspace, comment the part disp
        % and load parts 
        disp(['loading beads_analysis/prepro/beads_preproc_sub_', num2str(sub), '_block_', num2str(block)])
        load(['beads_analysis/prepro/beads_preproc_sub_', num2str(sub), '_block_', num2str(block)], 'data')
        
        first_data(block)  = data.firstdata;
        last_data(block)   = data.lastdata;
        
        
    end
    
    % Append data
    cfg             = [];
    all_firstdata   = ft_appenddata(cfg, first_data(1), first_data(2), first_data(3), first_data(4));
    all_lastdata    = ft_appenddata(cfg, last_data(1), last_data(2), last_data(3), last_data(4));
    
    % split data into all-but-last vs last easy and all-but-last vs last
    % diff
    thissub_first   = allsubs_allfirstERPs_keeptrials{1,sub};
    thissub_last    = allsubs_alllastERPs_keeptrials{1,sub};
    easy_first      = thissub_first{1,1};
    easy_last       = thissub_last{1,1};
    diff_first      = thissub_first{1,2};
    diff_last       = thissub_last{1,2};
    
    % first create the neighbours structure 
    cfg_neighb              = [];
    cfg_neighb.layout       = 'biosemi64.lay'; %in meters
    cfg_neighb.method       = 'distance';
    neighbours              = ft_prepare_neighbours(cfg_neighb, all_easydata);
    
    cfg                     = [];
    cfg.neighbours          = neighbours;   
                                            
    cfg.channel             = 'all';        
    cfg.latency             = [0 0.8];   
                                           

    % permutation tests
    cfg.method              = 'montecarlo';    
    cfg.statistic           = 'indepsamplesT'; 
                                               
    cfg.correctm            = 'cluster';
    cfg.clusteralpha        = 0.05;         
                                           
    cfg.clusterstatistic    = 'maxsum';                            
    cfg.minnbchan           = 2;            
    cfg.tail                = 0;           
    cfg.clustertail         = 0;
    cfg.alpha               = 0.025;        
    cfg.numrandomization    = 1000;       
    n_easyfirst             = size(easy_first.trial,1);
    n_easylast              = size(easy_last.trial,1);
    
    cfg.design              = [ones(1,n_easylast), ones(1,n_easyfirst)*2]; % design matrix
    cfg.ivar                = 1; % number or list with indices indicating the independent variable(s)

    % run stats 
    [stat]               = ft_timelockstatistics(cfg, easy_last, easy_first);
    
    % save stats 
    save(['beads_analysis/stats/beads_sub', num2str(sub), '_stats_1stlevel_easyfirst_vs_easylast'], 'stat')
    
    
    %%% -------------------------------------------------------------- %%%
    % Make a vector of all p-values associated with the clusters from ft_timelockstatistics.
    pos_cluster_pvals = [stat.posclusters(:).prob];

    % Then, find which clusters are deemed interesting to visualize, here we use a cutoff criterion based on the
    % cluster-associated p-value, and take a 5% two-sided cutoff (i.e. 0.025 for the positive and negative clusters,
    % respectively
    pos_clust = find(pos_cluster_pvals < 0.025);
    pos       = ismember(stat.posclusterslabelmat, pos_clust);


    % and now for the negative clusters...
    neg_cluster_pvals = [stat.negclusters(:).prob];
    neg_clust         = find(neg_cluster_pvals < 0.025);
    neg               = ismember(stat.negclusterslabelmat, neg_clust);

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
            print(gcf, '-dpng', ['beads_analysis/figures/stats_erps/beads_fig2_STAT_pos', num2str(i), '_subject', num2str(sub)])
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
            print(gcf, '-dpng', ['beads_analysis/figures/stats_erps/beads_fig2_STAT_neg', num2str(i), '_subject', num2str(sub)])
        end
    end
    %%% -------------------------------------------------------------- %%%
    
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
        cfg.saveaspng           = ['beads_analysis/figures/stats_erps/sub', num2str(sub), '_1stlevel_stats_easylast_vs_easyfirst'];
        cfg.toi                 = 0.2:0.01:0.8;
        cfg.layout              = 'biosemi64.lay';
        ft_clusterplot(cfg,stat); 
    end

    clear stat cfg cfg_neighb neigbours
    
    %% run stats for diff all-but-last vs last draws 
    
    % first create the neighbours structure 
    cfg_neighb              = [];
    cfg_neighb.layout       = 'biosemi64.lay'; %in meters
    cfg_neighb.method       = 'distance';
    neighbours              = ft_prepare_neighbours(cfg_neighb, all_diffdata);
    
    cfg                     = [];
    cfg.neighbours          = neighbours;   
                                            
    cfg.channel             = 'all';        
    cfg.latency             = [0 0.8];   
                                           

    % permutation tests
    cfg.method              = 'montecarlo';    
    cfg.statistic           = 'indepsamplesT'; 
                                               
    cfg.correctm            = 'cluster';
    cfg.clusteralpha        = 0.05;         
                                           
    cfg.clusterstatistic    = 'maxsum';                            
    cfg.minnbchan           = 2;            
    cfg.tail                = 0;           
    cfg.clustertail         = 0;
    cfg.alpha               = 0.025;        
    cfg.numrandomization    = 1000;       
    n_difffirst             = size(diff_first.trial,1);
    n_difflast              = size(diff_last.trial,1);
    
    cfg.design              = [ones(1,n_difflast), ones(1,n_difffirst)*2]; % design matrix
    cfg.ivar                = 1; % number or list with indices indicating the independent variable(s)

    % run stats 
    [stat]               = ft_timelockstatistics(cfg, diff_last, diff_first);
    
    % save stats 
    save(['beads_analysis/stats/beads_sub', num2str(sub), '_stats_1stlevel_difffirst_vs_difflast'], 'stat')
    
    %%% -------------------------------------------------------------- %%%
    % Make a vector of all p-values associated with the clusters from ft_timelockstatistics.
    pos_cluster_pvals = [stat.posclusters(:).prob];

    % Then, find which clusters are deemed interesting to visualize, here we use a cutoff criterion based on the
    % cluster-associated p-value, and take a 5% two-sided cutoff (i.e. 0.025 for the positive and negative clusters,
    % respectively
    pos_clust = find(pos_cluster_pvals < 0.025);
    pos       = ismember(stat.posclusterslabelmat, pos_clust);


    % and now for the negative clusters...
    neg_cluster_pvals = [stat.negclusters(:).prob];
    neg_clust         = find(neg_cluster_pvals < 0.025);
    neg               = ismember(stat.negclusterslabelmat, neg_clust);

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
            print(gcf, '-dpng', ['beads_analysis/figures/stats_erps/beads_fig3_STAT_pos', num2str(i), '_subject', num2str(sub)])
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
            print(gcf, '-dpng', ['beads_analysis/figures/stats_erps/beads_fig3_STAT_neg', num2str(i), '_subject', num2str(sub)])
        end
    end
    %%% -------------------------------------------------------------- %%%
    
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
        cfg.saveaspng           = ['beads_analysis/figures/stats_erps/sub', num2str(sub), '_1stlevel_stats_difflast_vs_difffirst'];
        cfg.toi                 = 0.2:0.01:0.8;
        cfg.layout              = 'biosemi64.lay';
        ft_clusterplot(cfg,stat); 
    end

    clear stat cfg cfg_neighb neigbours
    
    
end % end of subject loop

%% Within-Subjects statitsics ERPs 

% load the averaged ERPs 
% load('beads_analysis/erps/beads_allsubs_keeptrialsAllERPs', 'allsubERPallcond')
load('beads_analysis/erps/beads_allsubsTwoCondsERPs', 'allsubERPtwocond')
load('beads_analysis/erps/beads_allsubsAllButLastERPs', 'allsubfirstERPs')
load('beads_analysis/erps/beads_allsubsLastERPs', 'allsublastERPs')

%% Run group stats to diff vs easy trials 

% use the original data of one sub as example dataset to feed the cfg.neighbours matrix
% (to get sensor labels)
tmp_sub     = 1;
blocks      = 4;

% load the original easy and difficult data for tmp subject
for blocki = 1:blocks

    % load structure 
    disp(['loading beads_analysis/prepro/beads_preproc_sub_', num2str(tmp_sub), '_block_', num2str(blocki)])
    load(['beads_analysis/prepro/beads_preproc_sub_', num2str(tmp_sub), '_block_', num2str(blocki)], 'data')

    easy_data(blocki)   = data.easydata;
    diff_data(blocki)   = data.diffdata;

end

cfg                     = [];
all_easydata            = ft_appenddata(cfg, easy_data(1), easy_data(2), easy_data(3), easy_data(4));
all_diffdata            = ft_appenddata(cfg, diff_data(1), diff_data(2), diff_data(3), diff_data(4));

% split the averaged erps to diff and easy
subs = length(allsubERPtwocond);

for i = 1:subs
    
    avg_allsubs_easy{1,i} = allsubERPtwocond{1,i}{1,1};
    avg_allsubs_diff{1,i} = allsubERPtwocond{1,i}{1,2};
       
end

% first create the neighbours structure 
cfg_neighb              = [];
cfg_neighb.layout       = 'biosemi64.lay'; %in meters
cfg_neighb.method       = 'distance';
neighbours              = ft_prepare_neighbours(cfg_neighb, all_diffdata);

cfg                     = [];
cfg.neighbours          = neighbours;   

cfg.channel             = 'all';        
cfg.latency             = [0.3 0.6];   


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
cfg.numrandomization    = 'all';       

% create the design matrix
design                  = zeros(2, subs*2);
design(1,:)             = [1:subs 1:subs];
design(2,:)             = [ones(1,subs) ones(1,subs)*2];

cfg.design              = design;
cfg.uvar                = 1; % row 1 of design matrix contain var 1
cfg.ivar                = 2; % row 2 of design matrix contain var 2

[stat] = ft_timelockstatistics(cfg, avg_allsubs_diff{:}, avg_allsubs_easy{:});

% save stats 
save('beads_analysis/erps_group_stats/beads_stats_group_diff_vs_easy', 'stat')

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
        print(gcf, '-dpng', ['beads_analysis/figures/group_stats_erps/fig1_STAT_pos', num2str(i)])
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
        print(gcf, '-dpng', ['beads_analysis/figures/group_stats_erps/fig2_STAT_neg', num2str(i)])
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
    cfg.saveaspng           = 'beads_analysis/figures/group_stats_erps/beads_group_difflast_vs_difffirst';
    cfg.toi                 = 0.3:0.01:0.6; % times of interest
    cfg.layout              = 'biosemi64.lay';
    ft_clusterplot(cfg,stat); 
end

clear stat cfg cfg_neighb neigbours

%% Run Within-Subjects stats for all-but-last vs last darws easy trials

% split the averaged all-but-lust and last erps to diff and easy
subs = length(allsubERPtwocond);

for i = 1:subs
    
    avg_allsubs_first_easy{1,i} = allsubfirstERPs{1,i}{1,1};
    avg_allsubs_first_diff{1,i} = allsubfirstERPs{1,i}{1,2};
    
    avg_allsubs_last_easy{1,i} = allsublastERPs{1,i}{1,1};
    avg_allsubs_last_diff{1,i} = allsublastERPs{1,i}{1,2};
  
end

% first create the neighbours structure 
cfg_neighb              = [];
cfg_neighb.layout       = 'biosemi64.lay'; %in meters
cfg_neighb.method       = 'distance';
neighbours              = ft_prepare_neighbours(cfg_neighb, all_easydata);

cfg                     = [];
cfg.neighbours          = neighbours;   

cfg.channel             = 'all';        
cfg.latency             = [0 0.8];   


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
cfg.numrandomization    = 'all';       

% create the design matrix
design                  = zeros(2, subs*2);
design(1,:)             = [1:subs 1:subs];
design(2,:)             = [ones(1,subs) ones(1,subs)*2];

cfg.design              = design;
cfg.uvar                = 1; % row 1 of design matrix contain var 1
cfg.ivar                = 2; % row 2 of design matrix contain var 2

[stat] = ft_timelockstatistics(cfg, avg_allsubs_last_easy{:}, avg_allsubs_first_easy{:});

% save stats 
save('beads_analysis/erps_group_stats/beads_stats_group_easyfirst_vs_easylast', 'stat')

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
        print(gcf, '-dpng', ['beads_analysis/figures/group_stats_erps/fig3_STAT_pos', num2str(i)])
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
        print(gcf, '-dpng', ['beads_analysis/figures/group_stats_erps/fig4_STAT_neg', num2str(i)])
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
    cfg.saveaspng           = 'beads_analysis/figures/group_stats_erps/beads_group_easylast_vs_easyfirst';
    cfg.toi                 = 0.2:0.01:0.8;
    cfg.layout              = 'biosemi64.lay';
    ft_clusterplot(cfg,stat); 
end

clear stat cfg cfg_neighb neigbours

%% Run Within-Subjects stats for all-but-last vs last darws easy trials

% first create the neighbours structure 
cfg_neighb              = [];
cfg_neighb.layout       = 'biosemi64.lay'; %in meters
cfg_neighb.method       = 'distance';
neighbours              = ft_prepare_neighbours(cfg_neighb, all_diffdata);

cfg                     = [];
cfg.neighbours          = neighbours;   

cfg.channel             = 'all';        
cfg.latency             = [0 0.8];   


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
cfg.numrandomization    = 'all';       

% create the design matrix
design                  = zeros(2, subs*2);
design(1,:)             = [1:subs 1:subs];
design(2,:)             = [ones(1,subs) ones(1,subs)*2];

cfg.design              = design;
cfg.uvar                = 1; % row 1 of design matrix contain var 1
cfg.ivar                = 2; % row 2 of design matrix contain var 2

[stat] = ft_timelockstatistics(cfg, avg_allsubs_last_diff{:}, avg_allsubs_first_diff{:});

% save stats 
save('beads_analysis/erps_group_stats/beads_stats_group_difffirst_vs_difflast', 'stat')

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
        print(gcf, '-dpng', ['beads_analysis/figures/group_stats_erps/fig4_STAT_pos', num2str(i)])
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
        print(gcf, '-dpng', ['beads_analysis/figures/group_stats_erps/fig5_STAT_neg', num2str(i)])
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
    cfg.saveaspng           = 'beads_analysis/figures/group_stats_erps/beads_group_difflast_vs_difffirst';
    cfg.toi                 = 0.2:0.01:0.8;
    cfg.layout              = 'biosemi64.lay';
    ft_clusterplot(cfg,stat); 
end

clear stat cfg cfg_neighb neigbours

%% Run TFR statistics (between-trials/first level)

% load the TFR data

load('beads_analysis/tfr/beads_allsubs_tfrEasy', 'allsub_TFReasy')
load('beads_analysis/tfr/beads_allsubs_tfrDiff', 'allsub_TFRdiff')
load('beads_analysis/tfr/beads_allsubs_tfrFirst', 'allsub_TFRfirst')
load('beads_analysis/tfr/beads_allsubs_tfrLast', 'allsub_TFRlast')

% Run statitics at the trial-level for each subject first 
% Cluster-based permutation tests will be used to look at differences
% between conditions: (easy vs diff) and draws: (all-but-last vs last
% darw)with a focus at beta band activity from TFR analysis

% extract current subject mat file from a given cell loaded above and run
% analysis
for sub = 1:nsubs 
    
    % load the original easy and difficult data for current sub 
    for blocki = 1:blocks
        
        % load structure 
        disp(['loading beads_analysis/prepro/beads_preproc_sub_', num2str(sub), '_block_', num2str(blocki)])
        load(['beads_analysis/prepro/beads_preproc_sub_', num2str(sub), '_block_', num2str(blocki)], 'data')

        easy_data(blocki)   = data.easydata;
        diff_data(blocki)   = data.diffdata;

    end
    
    cfg                     = [];
    all_easydata            = ft_appenddata(cfg, easy_data(1), easy_data(2), easy_data(3), easy_data(4));
    all_diffdata            = ft_appenddata(cfg, diff_data(1), diff_data(2), diff_data(3), diff_data(4));
    
    % run cluster-based permutations using the monte-carlo method for easy
    % vs difficult trials/epochs for each subject seperately 
    thissub_TFReasy         = allsub_TFReasy{1,sub};
    thissub_TFRdiff         = allsub_TFRdiff{1,sub};
    
    % permutation tests for TFR data run the same way as for ERP data
    
    % create the neighbours structure which will be critical to compute
    % stats
    cfg_neighb              = [];
    cfg_neighb.layout       = 'biosemi64.lay'; %in meters
    cfg_neighb.method       = 'distance';
    neighbours              = ft_prepare_neighbours(cfg_neighb, all_easydata);
    
    cfg                     = [];
    cfg.neighbours          = neighbours;   % the neighbours specify for each sensor with
                                            % which other sensors it can form clusters
    cfg.channel             = 'all';        % cell-array with selected channel labels
    cfg.latency             = 'all';        % time interval over which the experimental - I will use 'all' for now
                                            % conditions must be compared (in seconds)
                                            
    cfg.frequency           = 30;           % freq of interest 
    
    % permutation tests
    cfg.method              = 'montecarlo';     % use the Monte Carlo Method to calculate the significance probability
    cfg.statistic           = 'ft_statfun_indepsamplesT';  % use the independent samples T-statistic as a measure to
                                                % evaluate the effect at the sample level
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

    % design matrix
    design                  = zeros(1,size(thissub_TFRdiff.powspctrm,1) + size(thissub_TFReasy.powspctrm,1));
    
    design(1,1:size(thissub_TFRdiff.powspctrm,1)) = 1;
    design(1,(size(thissub_TFRdiff.powspctrm,1)+1):(size(thissub_TFRdiff.powspctrm,1)+size(thissub_TFReasy.powspctrm,1))) = 2;
    
    cfg.design              = design; 
    cfg.ivar                = 1;
    
    [stat]                  = ft_freqstatistics(cfg, thissub_TFRdiff, thissub_TFReasy);
    
 
end

%% Run TFR statistics (Within-Subjects/second level)

% load Grand Average TFRs with individual trials (the power spectrum structure must be 4-D)
load('beads_analysis/tfr/beads_ga_TFReasy_ind', 'ga_TFReasy_ind')
load('beads_analysis/tfr/beads_ga_TFRdiff_ind', 'ga_TFRdiff_ind')
load('beads_analysis/tfr/beads_ga_TFRfirst_ind', 'ga_TFRfirst_ind')
load('beads_analysis/tfr/beads_ga_TFRlast_ind', 'ga_TFRlast_ind')
load('beads_analysis/tfr/beads_ga_TFRfirst_easy_ind', 'ga_TFRfirst_easy_ind')
load('beads_analysis/tfr/beads_ga_TFRlast_easy_ind', 'ga_TFRlast_easy_ind')
load('beads_analysis/tfr/beads_ga_TFRfirst_diff_ind', 'ga_TFRfirst_diff_ind')
load('beads_analysis/tfr/beads_ga_TFRlast_diff_ind', 'ga_TFRlast_diff_ind')


%% Run TFR statistics for difficult vs easy trials

% create the neighbours structure which will be critical to compute stats
cfg_neighb              = [];
cfg_neighb.layout       = 'biosemi64.lay'; %in meters
cfg_neighb.method       = 'distance';
neighbours              = ft_prepare_neighbours(cfg_neighb, ga_TFRdiff_ind);

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
subj                    = 8; % this needs to be changed everytime we add a subject
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

[stat] = ft_freqstatistics(cfg, ga_TFRdiff_ind, ga_TFReasy_ind);

% save stats 
save('beads_analysis/tfr_group_stats/beads_TFRstats_group_diff_vs_easy', 'stat')

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
    cfg.saveaspng           = 'beads_analysis/figures/group_stats_tfrs/beads_groupTFR_diff_vs_easy';
    %cfg.toi                 = 0.0:0.01:0.6; % times of interest
    cfg.layout              = 'biosemi64.lay';
    ft_clusterplot(cfg,stat); 
end

clear stat cfg cfg_neighb

%% Run TFR statistics for easy last vs all-but-last draws

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
subj                    = 8; % this needs to be changed everytime we add a subject
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

[stat] = ft_freqstatistics(cfg, ga_TFRlast_easy_ind, ga_TFRfirst_easy_ind);

% save stats 
save('beads_analysis/tfr_group_stats/beads_TFRstats_group_last_vs_first_easy', 'stat')

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
    cfg.saveaspng           = 'beads_analysis/figures/group_stats_tfrs/beads_groupTFR_last_vs_first_easy';
    %cfg.toi                 = 0.0:0.01:0.6; % times of interest
    cfg.layout              = 'biosemi64.lay';
    ft_clusterplot(cfg,stat); 
end

clear stat cfg cfg_neighb

%% Run TFR statistics for difficult last vs all-but-last draws

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
subj                    = 8; % this needs to be changed everytime we add a subject
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

[stat] = ft_freqstatistics(cfg, ga_TFRlast_diff_ind, ga_TFRfirst_diff_ind);

% save stats 
save('beads_analysis/tfr_group_stats/beads_TFRstats_group_last_vs_first_diff', 'stat')

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
    cfg.saveaspng           = 'beads_analysis/figures/group_stats_tfrs/beads_groupTFR_last_vs_first_diff';
    %cfg.toi                 = 0.0:0.01:0.6; % times of interest
    cfg.layout              = 'biosemi64.lay';
    ft_clusterplot(cfg,stat); 
end

clear stat cfg cfg_neighb




