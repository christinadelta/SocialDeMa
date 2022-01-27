%% BEADS FORMAL PREPROCESSING AND ANALYSES SCRIPT W/ FIELDTRIP

% christinadelta 
% Date: Jan 2022

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

%% Define paths and variables 

% Define paths 
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
nconds              = 4;
conditions          = {'blueEasy', 'greenEasy', 'blueDiff', 'greenDiff'};

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
        
        % re-write the trl list to the cfg struct 
        cfg.trl                     = trl;

        clear trialend trialstart j i counter cnt trl trialnum 
        
        %% re-reference/preprocess

        cfg.reref           = 'yes';
        cfg.refchannel      = {'EXG1' 'EXG2'};
        cfg.demean          = 'yes';
        cfg.baselinewindow  = [-0.2 0];
        data                = ft_preprocessing(cfg);
        
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

    % Load the matfiles and concatenate all blocks 
    for block = 1:blocks
        
        disp(['loading beads_analysis/prepro/beads_preproc_sub_', num2str(subI), '_block_', num2str(block)])
        load(['beads_analysis/prepro/beads_preproc_sub_', num2str(subI), '_block_', num2str(block)], 'data')
        
        partdata(block) = data;
        
    end
    
    % Append data
    cfg     = [];
    data    = ft_appenddata(cfg, partdata(1), partdata(2), partdata(3), partdata(4));
    
    %% Remove artifacts with ICA 
    
    % THIS PART IS PROBABLY NOT NEEDED 
    % COMMENT IT IF ICA WONT RUN 
    
    % In general, in this experiment we don't need to run ICA for blink
    % removal because we instruct participants NOT to blink during this
    % time. Possibly will only need to filter.
    
    % Detect eog artifacts using ICA
    cfg             = [];
    cfg.method      = 'runica';
    cfg.channel     = 1:64; % EEG channels only
    datacomp        = ft_componentanalysis(cfg, data);
    % save('beads_analysis/ica/datacomp', 'datacomp')

    % plot the components to detect the artifacts
    figure
    k = 1; f = 1;

    for icomp=1:length(datacomp.topo)

      if k>20
        k = 1;
        figure
      end

      cfg           = [];
      cfg.layout    = 'biosemi64.lay';
      cfg.xlim      = [icomp icomp];

      subplot(4,5,k);
      ft_topoplotER(cfg, datacomp);
      title(icomp);

      k = k+1;
    end

    % remove components that reflect eog artifacts
    cfg             = [];
    cfg.component   = [21 42 63 68]; % the exact numbers varies per run
    clean_data      = ft_rejectcomponent(cfg, datacomp);
    
    %% filter data (clean) and run timelock (ERPs) analysis (3 analyses)
    
    % 1) Run ERPs analysis: 1 average ERP for each sequence
    % Average ERP: 1:end (epochs) ( 52 averaged ERPs total)  
    for iBlock = 1:blocks 
        for itrial = 1:blocktrials
            cfg                         = [];
            cfg.lpfilter                = 'yes';
            cfg.lpfreq                  = 40;
            cfg.trials                  = find(data.trialinfo(:,3) == iBlock & data.trialinfo(:,2) == itrial);
            timelock{iBlock, itrial}    = ft_timelockanalysis(cfg, data);
            
            % baseline correction
            cfg                         = [];
            cfg.baseline                = [-0.2 0];
            timelock{iBlock, itrial}    = ft_timelockbaseline(cfg, timelock{iBlock, itrial});
        end
    end
    
    % save the erp analysis 
    
    % plot the ERPs over all sensors, over parietal sensors, over frontal
    % sensors 
    
    
    
    % 2) Run ERPs analysis: For every condition, in each "sequence" all epochs/draws are a separate ERP 
    % [e.g., if block 1 trial 1 has 4 draws/epochs of cond 1 -- we compute 4 seperate ERPs]  
    for cond = 1:nconds
        for itrial = 1:totaltrials
            cfg                         = [];
            cfg.lpfilter                = 'yes';
            cfg.lpfreq                  = 40;
            cfg.trials                  = find(data.trialinfo(:,1) == cond & data.trialinfo(:,2) == itrial);
            timelock{cond, itrial}    = ft_timelockanalysis(cfg, data);

            % baseline correction
            cfg                         = [];
            cfg.baseline                = [-0.2 0];
            timelock{cond, itrial}    = ft_timelockbaseline(cfg, timelock{iBlock, itrial});
        end
    end
    
    
    
    
    
    
    
    
    
    
    
    % 3) Run ERPs analysis: 2 average ERPs for each sequence/condition 
    % ERP 1: 1:end-1 & ERP 2: end (last epoch)
    for iBlock = 1:blocks 
        for itrial = 1:blocktrials
            for cond = 1:nconds 

                cfg = [];
                cfg.lpfilter = 'yes';
                cfg.lpfreq = 40;
                cfg.trials = find(data.trialinfo(:,3) == iBlock & data.trialinfo(:,2) == itrial & data.trialinfo(:,1) == cond);
                if ~isempty(cfg.trials)
                    timelock{iBlock, itrial, cond} = ft_timelockanalysis(cfg, data)
                end
            end
        end
    end
    
   
 
    %% 
    
    
    %% Compute statistics (contrasts)
    
    
  
end % end of subject loop






