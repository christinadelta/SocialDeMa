%% BEADS FORMAL PREPROCESSING AND ANALYSES SCRIPT W/ FIELDTRIP

% christinadelta 
% Date: Jan 2022
% last Update 10/02/2022

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
% 1). Visualise ERPs in different ways (check)
% 2). Run stats
% 3). Extract Frequences 

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
        
        cfg.trl             = trl_easy;
        easy_data           = ft_preprocessing(cfg);
        redata.easy_data    = easy_data;
        
        cfg.trl             = trl_diff;
        diff_data           = ft_preprocessing(cfg);
        redata.diff_data    = diff_data;
        
        cfg.trl             = first_trl;
        first_data          = ft_preprocessing(cfg);
        redata.first_data   = first_data;
        
        cfg.trl             = last_trl;
        last_data           = ft_preprocessing(cfg);
        redata.last_data    = last_data;
        
        % only keep eeg channels from now on % for all data structures 
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
        
        disp(['loading beads_analysis/prepro/beads_preproc_sub_', num2str(subI), '_block_', num2str(block)])
        load(['beads_analysis/prepro/beads_preproc_sub_', num2str(subI), '_block_', num2str(block)], 'data')
        
        partdata(block) = data.alldata;
        
        % visualise the block trials/epochs
%         cfg                 = [];
%         ft_databrowser(cfg, partdata(block))
        
    end
    
    % Append data
    cfg     = [];
    alldata = ft_appenddata(cfg, partdata(1), partdata(2), partdata(3), partdata(4));
    
    
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
    
    % I will run 4 different ERP analyses only to be able to visualise them
    % properly and to explore the data. 
    
    % The ERPs that will mainly be used from now on are: [averagedERPsSecond, timelock_all]
    
    % 1) Run ERPs analysis: 1 average ERP for each sequence
    % Average ERP: 1:end (epochs) ( 52 averaged ERPs total) - 4 Conditions   
    for iCond = 1:nconds
        for iBlock = 1:blocks
            for iTrial = 1:blocktrials
                cfg                                 = [];
                cfg.lpfilter                        = 'yes';
                cfg.lpfreq                          = 40;
                tmp                                 = find(alldata.trialinfo(:,1) == iCond & alldata.trialinfo(:,3) == iBlock & alldata.trialinfo(:,2) == iTrial);
                
                if ~isempty(tmp) % if tmp is not empty run timelock analysis on current trial/sequence
                    cfg.trials                      = tmp;
                    timelock{iCond,iBlock,iTrial}   = ft_timelockanalysis(cfg, alldata);
                    
                    % baseline correction
                    cfg                             = [];
                    cfg.baseline                    = [-0.2 0];
                    timelock{iCond,iBlock,iTrial}   = ft_timelockbaseline(cfg, timelock{iCond,iBlock,iTrial});
                    
                end
            end % end of trials loop
            
        end % end of blocks loop
        
        % remove empty cells from timelock cell
        timelock(cellfun(@(timelock) any(isempty(timelock)),timelock)) = [];
        
        % save timelock for each condition in new cell and clear it (otherwise timelock
        % cell is flattened using the removal function and I can't separate
        % the conditions)
        averagedERPs{iCond,:} = timelock; clear timelock tmp 
        
    end % end of conds loop
    
    % 2) Run ERPs analysis: 1 average ERP for each sequence
    % Average ERP: 1:end (epochs) ( 52 averaged ERPs total) - 2 Conditions
    totalconds  = 2; 
    
    % add a column of ones and twos for easy difficukt trials in
    % data.trialinfo
    trlinfo         = alldata.trialinfo;
    totaldraws      = length(trlinfo);
    condtrials      = blocktrials * totalconds;
    
    for i = 1:totaldraws
        
        if trlinfo(i,1) == 1 | trlinfo(i,1) == 2
            
            trlinfo(i,4) = 1; % add easy index 
            
        elseif trlinfo(i,1) == 3 | trlinfo(i,1) == 4
            
            trlinfo(i,4) = 2; % add difficult index 
        end
    end
    
    data.trialinfo = trlinfo; clear tlrinfo 
    
    % run third ERP analysis using 2 conditions (easy vs diff for each
    % sequence)
    for thiscond = 1:totalconds
        for block = 1:blocks
            for trial = 1:blocktrials
                
                cfg                                 = [];
                cfg.lpfilter                        = 'yes';
                cfg.lpfreq                          = 40;
                tmp                                 = find(alldata.trialinfo(:,4) == thiscond & alldata.trialinfo(:,3) == block & alldata.trialinfo(:,2) == trial);
                
                if ~isempty(tmp)
                    
                    cfg.trials                      = tmp;
                    timelock{thiscond,block,trial}  = ft_timelockanalysis(cfg, alldata);
                    
                    % baseline correction
                    cfg                             = [];
                    cfg.baseline                    = [-0.2 0];
                    timelock{thiscond,block,trial}  = ft_timelockbaseline(cfg, timelock{thiscond,block,trial});
  
                end % end of if
            end % end of trials
        end % end of blocks
        
        % remove empty cells from timelock cell
        timelock(cellfun(@(timelock) any(isempty(timelock)),timelock)) = [];
        
        % save timelock for each condition in new cell and clear it (otherwise timelock
        % cell is flattened using the removal function and I can't separate
        % the conditions)
        averagedERPSecond{thiscond,:} = timelock; clear timelock tmp % this is the second averaged ERP analysis whith 2 conds
        
    end % end of conds loop
    
    % 3) Run 3rd ERP analysis where you average all diff epochs and all easy
    % epochs. This will probs be used only for visualising (2 epochs/trials- easy/diff). 
    for thiscond = 1:totalconds 
        
        cond_all    = 0;
        cnt         = 0;
        
        for block = 1:blocks
            for trial = 1:blocktrials
                
                cfg                                 = [];
                cfg.lpfilter                        = 'yes';
                cfg.lpfreq                          = 40;
                tmp                                 = find(alldata.trialinfo(:,4) == thiscond & alldata.trialinfo(:,3) == block & alldata.trialinfo(:,2) == trial);
                
                if ~isempty(tmp)
                    
                    % store tmp indexes of this cond/trial to cond_all
                    tmp_length                      = length(tmp);
                    cond_all(cnt+1:cnt+tmp_length)  = tmp;
                    cnt                             = cnt + tmp_length;
                end 
            end % end of trials loop
        end % end of blocks loop
        
        cfg.trials                  = cond_all;
        timelock{thiscond}              = ft_timelockanalysis(cfg, alldata);
        clear cond_all cnt tmp_length tmp
    
    end % end of conds loop
    
    twocondsERPs_avr                = timelock; clear timelock

    % 4) Run ERPs analysis: 1 averaged ERP for draws 1:end-1 (first draw until last-1) & 1 averaged
    % ERP end (last draw) across all sequences/trials for each condition
    % (easy, difficult). For this we first need to make
    % sure that for each condition (easy, diff) we get 2 ERPs; so, let's
    % first split the draws in different matrices. One that will contain
    % all the [1: end-1] draws and one matrix that will contain all the
    % [end] draws. Then run timelock analysis (ERPs) and store in cell for
    % each condition separately.
    
    for thiscond = 1:totalconds
        
        % init variables
        tmp_first                       = 0;
        tmp_last                        = 0;
        c                               = 0; % counter index
        l                               = 1; % last draw index
    
        for block = 1:blocks
            for trial = 1:blocktrials
                
                tmp                     = find(alldata.trialinfo(:,4) == thiscond & alldata.trialinfo(:,3) == block & alldata.trialinfo(:,2) == trial);
                
                % if tmp is not empty for this cond/block/trial, split the
                % draws/epochs in [1:end-1] and [end]
                if ~isempty(tmp)
                    
                    t                   = length(tmp)-1; 
                    tmp_first(c+1:c+t)  = tmp(1:end-1); 
                    tmp_last(:,l)       = tmp(end);
                    
                    % update counter and l
                    c                   = c + t;
                    l                   = l + 1;
                    
                end
            end % end of trials loop
        end % end of blocks loop
        
        % run timelock/ERP analysis on the trials [1:end-1] for
        % this condition and this block
        cfg.trials                      = tmp_first;
        timelock_first{thiscond}        = ft_timelockanalysis(cfg, alldata);

        % baseline correction
        cfg                             = [];
        cfg.baseline                    = [-0.2 0];
        timelock_first{thiscond}        = ft_timelockbaseline(cfg, timelock_first{thiscond});

        % run timelock/ERP analysis on the trials [end/last] for
        % this condition and this block
        cfg.trials                      = tmp_last;
        timelock_last{thiscond}         = ft_timelockanalysis(cfg, alldata);

        % baseline correction
        cfg                             = [];
        cfg.baseline                    = [-0.2 0];
        timelock_last{thiscond}         = ft_timelockbaseline(cfg, timelock_last{thiscond});
        
        % clear temporal variables to re-initialise them for the next
        % condition 
        clear tmp_first tmp_last tmp c l t 

    end % end of conditions loop
    
    timelock_all{1,:}   = timelock_first;
    timelock_all{2,:}   = timelock_last;
    
    
    clear timelock_first timelock_last allButLast_one allButLast_two last_one last_two
    
    % save the averaged [1:end-1] and [last] erp analysis for each subject 
    save(['beads_analysis/erps/beads_sub_', num2str(subI), '_alldraws_erps'], 'timelock_all')
   
    % save erp analysis (4 conds)
    save(['beads_analysis/erps/beads_sub_', num2str(subI),  '_allerps'], 'averagedERPs')
    
    % save erp analysis (2 conds)
    save(['beads_analysis/erps/beads_sub_', num2str(subI),  '_allerps_second'], 'averagedERPSecond')
    
    % save averaged erps (2 conds) analysis
    save(['beads_analysis/erps/beads_sub_', num2str(subI),  '_erps_avr'], 'twocondsERPs_avr')
    
    %% Plot the ERPs 
    
    % load the data (if needed) and plot the ERPs over all sensors, over parietal sensors, over frontal
    % sensors 
    % load(['beads_analysis/erps/beads_sub_', num2str(subI),  '_erps_averaged_'], 'averagedERPs');
    
    % define new variables 
    erps        = {'AllButLast', 'Last'};
    
%    % plot [all but last and last draw] ERPs for each condition over all sensors
%     figure
%     for i = 1:2
%         
%         erp_tmp = timelock_all{i};
%         
%         subplot(2,2,i)
%         ft_singleplotER([], erp_tmp{1,1}, erp_tmp{1,2});
%         title(erps{i});
%         legend('easy', 'difficult')
%     end
    
    % ------------------------------------
    % a) plot over frontal sensors (diff vs easy) averaged ERPs
    cfg = [];
    cfg.channel = {'F1', 'F3', 'F5', 'F7', 'Fz', 'F2', 'F4', 'F6', 'F8'};
    figure; 
    ft_singleplotER(cfg, twocondsERPs_avr{1,1}, twocondsERPs_avr{1,2})
    title('forntal channels averaged ERPs')
    legend('easy', 'difficult')
    
    print(gcf, '-dpng', ['beads_analysis/figures/beads_sub', num2str(subI), '_fig1_avrERP'])

    % b) plot over parietal sensors (diff vs easy) averaged ERPs
    cfg = [];
    cfg.channel = {'P1', 'P3', 'P5', 'P7', 'P9', 'Pz', 'P2', 'P4', 'P6', 'P8', 'P10'};
    figure; 
    ft_singleplotER(cfg, twocondsERPs_avr{1,1}, twocondsERPs_avr{1,2})
    title('parietal channels averaged ERPs')
    legend('easy', 'difficult')
    
    print(gcf, '-dpng', ['beads_analysis/figures/beads_sub', num2str(subI), '_fig2_avrERP'])

    % -------------------------------------------
    % c) plot [all draws - last] vs last draw for easy cond (over frontal
    % electrodes)
    cfg = [];
    cfg.channel = {'F1', 'F3', 'F5', 'F7', 'Fz', 'F2', 'F4', 'F6', 'F8'};
    figure; 
    ft_singleplotER(cfg, timelock_all{1,1}{1,1}, timelock_all{2,1}{1,1})
    title('all draws but last vs last draw easy-frontal')
    legend('allButLast', 'last')
    
    % save figures
    print(gcf, '-dpng', ['beads_analysis/figures/beads_sub', num2str(subI), '_frontal_easy_mainERPs'])

    % d) plot [all draws - last] vs last draw for easy cond (over parietal
    % electrodes)
    cfg = [];
    cfg.channel = {'P1', 'P3', 'P5', 'P7', 'P9', 'Pz', 'P2', 'P4', 'P6', 'P8', 'P10'};
    figure; 
    ft_singleplotER(cfg, timelock_all{1,1}{1,1}, timelock_all{2,1}{1,1})
    title('all draws but last vs last draw easy-parietal')
    legend('allButLast', 'last')
    
    print(gcf, '-dpng', ['beads_analysis/figures/beads_sub', num2str(subI), '_parietal_easy_mainERPs'])
    
    % e) plot [all draws - last] vs last draw for diff cond (over frontal
    % electrodes)
    cfg = [];
    cfg.channel = {'F1', 'F3', 'F5', 'F7', 'Fz', 'F2', 'F4', 'F6', 'F8'};
    figure; 
    ft_singleplotER(cfg, timelock_all{1,1}{1,2}, timelock_all{2,1}{1,2})
    title('all draws but last vs last draw diff-frontal')
    legend('allButLast', 'last')
    
    print(gcf, '-dpng', ['beads_analysis/figures/beads_sub', num2str(subI), '_frontal_diff_mainERPs'])
    
    % f) plot [all draws - last] vs last draw for diff cond (over parietal
    % electrodes)
    cfg = [];
    cfg.channel = {'P1', 'P3', 'P5', 'P7', 'P9', 'Pz', 'P2', 'P4', 'P6', 'P8', 'P10'};
    figure; 
    ft_singleplotER(cfg, timelock_all{1,1}{1,2}, timelock_all{2,1}{1,2})
    title('all draws but last vs last draw diff-parietal')
    legend('allButLast', 'last')
    
    print(gcf, '-dpng', ['beads_analysis/figures/beads_sub', num2str(subI), '_parietal_diff_mainERPs'])
    
    % ----------------------------------------------------------
    % plot ERPs on interactive mode (plot all but last and last for easy cond)
    cfg = [];
    cfg.layout  = 'biosemi64.lay';
    cfg.interactive = 'yes';
    figure; ft_multiplotER(cfg, timelock_all{1,1}{1,1}, timelock_all{2,1}{1,1})
    
    % plot ERPs on interactive mode (plot all but last and last draws for diff cond)
    cfg = [];
    cfg.layout  = 'biosemi64.lay';
    cfg.interactive = 'yes';
    figure; ft_multiplotER(cfg, timelock_all{1,1}{1,2}, timelock_all{2,1}{1,2})
    
    % -----------------------------------------------------------------
    % maybe plot main ERPs in topoplots?
    % ALL but last, EASY
    cfg = [];
    cfg.xlim = [0.3 0.5];
    % cfg.zlim = [0 6e-14]; % this command messes up with the plot for some reason 
    cfg.layout  = 'biosemi64.lay';
    cfg.parameter = 'avg';
    cfg.interactive = 'yes';
    figure; ft_topoplotER(cfg, timelock_all{1,1}{1,1}); colorbar
    
    % save
    print(gcf, '-dpng', ['beads_analysis/figures/beads_sub', num2str(subI), '_allbutlast_easy_topo'])
    
    % Last, Easy
    cfg = [];
    cfg.xlim = [0.3 0.5];
    cfg.layout  = 'biosemi64.lay';
    cfg.parameter = 'avg';
    cfg.interactive = 'yes';
    figure; ft_topoplotER(cfg, timelock_all{2,1}{1,1}); colorbar
    
    print(gcf, '-dpng', ['beads_analysis/figures/beads_sub', num2str(subI), '_last_easy_topo'])
    
    % ALL but last, DIFF
    cfg = [];
    cfg.xlim = [0.3 0.5];
    cfg.layout  = 'biosemi64.lay';
    cfg.parameter = 'avg';
    cfg.interactive = 'yes';
    figure; ft_topoplotER(cfg, timelock_all{1,1}{1,2}); colorbar
    
    % save
    print(gcf, '-dpng', ['beads_analysis/figures/beads_sub', num2str(subI), '_allbutlast_diff_topo'])
    
    % Last, DIFF
    cfg = [];
    cfg.xlim = [0.3 0.5];
    cfg.layout  = 'biosemi64.lay';
    cfg.parameter = 'avg';
    cfg.interactive = 'yes';
    figure; ft_topoplotER(cfg, timelock_all{2,1}{1,2}); colorbar
    
    % save
    print(gcf, '-dpng', ['beads_analysis/figures/beads_sub', num2str(subI), '_last_diff_topo'])

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
        disp(['loading beads_analysis/prepro/beads_preproc_sub_', num2str(subI), '_block_', num2str(block)])
        load(['beads_analysis/prepro/beads_preproc_sub_', num2str(subI), '_block_', num2str(block)], 'data')

        easy_data(blocki)   = data.easydata;
        diff_data(blocki)   = data.diffdata;
        first_data(blocki)  = data.firstdata;
        last_data(blocki)   = data.lastdata;
        
    end
    
    cfg             = [];
    all_easydata    = ft_appenddata(cfg, easy_data(1), easy_data(2), easy_data(3), easy_data(4));
    all_diffdata    = ft_appenddata(cfg, diff_data(1), diff_data(2), diff_data(3), diff_data(4));
    all_firstdata   = ft_appenddata(cfg, first_data(1), first_data(2), first_data(3), first_data(4));
    all_lastdata    = ft_appenddata(cfg, last_data(1), last_data(2), last_data(3), last_data(4));
    
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
    TFRwave_easy   = ft_freqanalysis(cfg, all_easydata);
    TFRwave_diff   = ft_freqanalysis(cfg, all_diffdata);
%     
    
    cfg = [];
    cfg.baseline     = [-0.2 0];
    cfg.baselinetype = 'absolute';
    cfg.marker       = 'on';
    cfg.showlabels   = 'yes';
    cfg.layout       = 'biosemi64.lay';
    figure; ft_multiplotTFR(cfg, TFRwave_easy); ft_multiplotTFR(cfg, TFRwave_diff);
    
    % ---------------------------
    % run TFR for Group 2 conditions
    cfg = [];
    cfg.channel    = 'all';
    cfg.method     = 'wavelet';
    cfg.width      = 7; % width definition is important 
    cfg.pad        = 10;
    cfg.output     = 'pow';
    cfg.foi        = 1:1:30;
    cfg.toi        = 'all'; % for computation efficiency use 'all'
    TFRwave_first   = ft_freqanalysis(cfg, all_firstdata);
    TFRwave_last   = ft_freqanalysis(cfg, all_lastdata);
%     
    
    cfg = [];
    cfg.baseline     = [-0.2 0];
    cfg.baselinetype = 'absolute';
    cfg.marker       = 'on';
    cfg.showlabels   = 'yes';
    cfg.layout       = 'biosemi64.lay';
    figure; ft_multiplotTFR(cfg, TFRwave_first); ft_multiplotTFR(cfg, TFRwave_last);

   

   
    
  
    %% Compute contrasts 
    
    % first compute contrasts for averaged trials/sequences (easy vs diff)
    
    
    
    
    %% Compute statistics (contrasts)
    
    
  
end % end of subject loop

