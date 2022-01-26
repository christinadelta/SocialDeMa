%% ANALYSIS OF THE BEADS TASK

% analysis of the EEG data (beads task) using the FieldTrip toolbox

% TODO:
% 1. loop over blocks (every block is a seperate .bdf file)
% 2. loop over subjects

% define subject path 
sub_path = '/Users/christinadelta/Desktop/os_data/beads';
sub_data = fullfile(sub_path, 'sub_08_beads_block_01.bdf');

% check the biosemi64 layout that will be used for plotting topomaps
cfglayout               = [];
cfglayout.layout        = 'biosemi64.lay';
lo                      = ft_layoutplot(cfglayout);

%% load the data  and define trials 
cfg                     = [];
cfg.dataset             = sub_data;

cfg.trialdef.eventtype  = 'STATUS';
cfg.trialdef.eventvalue = [1 2 3 4 102 103];
cfg.trialdef.prestim    = 0.2;
cfg.trialdef.poststim   = 0.8;
cfg                     = ft_definetrial(cfg);

% extract trl list
trl                     = cfg.trl;

tstart                  = 102;
tend                    = 103;

blocktrials             = length(find(trl(:,4) == tstart));
trialstart              = find(trl(:,4) == tstart);
trialend                = find(trl(:,4) == tend);

counter                 = 0;

% add trial number 
for i = 1:blocktrials
    
    tmp                 = length(trialstart(i)+1: trialend(i)-1);
    for j = 1:tmp
        
        cnt             = counter + j;
        trialnum(cnt,1) = i;
    end
    
    counter             = counter + tmp;
    clear tmp
  
end

% remove trialstart and trialend from the list
trl(trl(:,4) == tstart, :)  = [];
trl(trl(:,4) == tend, :)    = [];
trl(:,5)                    = trialnum; % add trialnum to the main list

% re-write the trl list to the cfg struct 
cfg.trl                     = trl;

clear trialend trialstart j i counter cnt trl trialnum 

%% re-reference/preprocess

cfg.reref           = 'yes';
cfg.refchannel      = {'EXG1' 'EXG2'};
cfg.demean          = 'yes';
cfg.baselinewindow  = [-0.2 0];
data                = ft_preprocessing(cfg);

% save the preprocessed file for now 
save(['beads_analysis/prepro/beads_prepro'], 'data')

%% visualise the trials/epochs of interest 

% NOTE: in typical analysis, each trial includes one stimulus and we epoch
% the EEG signal around this stimulus. However, in the beads task each
% trial contains a number of stimuli/draws (up to 10), thus here we define trials
% and we epoch the EEG signal around each stimulus/draw in that sequnce.
% Thus, we have 13 trials/sequences in each block, and length(trl) epochs. 

% we can either visualise the epochs for each trial/sequence seperately or
% all together. For simplicity, I'll visualise them together. 
cfg                 = [];
ft_databrowser(cfg, data)

% maybe also visualise the epochs per trial/sequence?

%% artifact rejection 

% 1. detect eog artifacts using ICA
cfg             = [];
cfg.method      = 'runica';
cfg.channel     = 1:64; % EEG channels only
datacomp        = ft_componentanalysis(cfg, data);
save('beads_analysis/ica/datacomp', 'datacomp')

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

save('beads_analysis/data_clean', 'clean_data') % save the clean data

% 2. remove artifacts manually 
% remove the mean
cfg             = [];
cfg.demean      = 'yes';
data            = ft_preprocessing(cfg, data);

% shuffle the trials for non-biased artifact rejection
shuffle         = randperm(length(data.trial));
datashuff       = data;

for i = 1:length(data.trial)
    
  datashuff.trial{i}        = data.trial{shuffle(i)};
  datashuff.time{i}         = data.time{shuffle(i)};
  datashuff.trialinfo(i,:)  = data.trialinfo(shuffle(i),:);
  
end

% browse for artifact rejection
cfg             = [];
cfg.channel     = 'EEG';
cfg.continuous  = 'no';
cfg             = ft_databrowser(cfg, clean_data);
data            = ft_rejectartifact(cfg,data);

% visual artifact rejection in summary mode
cfg             = [];
cfg.method      = 'summary';
data            = ft_rejectvisual(cfg, data);

% save('analysis/data_clean2', 'data')

%% timelock analysis 

% define conditions 
conds           = {'easy_blue', 'easy_green', 'difficult_blue', 'difficult_green'};
conditions      = length(conds);
% how many trials/sequences per block?
blocktrials     = 13;

% load the clean data
% load('beads_analysis/data_clean', 'clean_data')

% compute timelocked averages for each trial
for i = 1:blocktrials

    cfg             = [];
    cfg.lpfilter    = 'yes';
    cfg.lpfreq      = 40;
    cfg.trials      = find(clean_data.trialinfo(:,2) == i);
    timelock{i}     = ft_timelockanalysis(cfg, clean_data);
    
    % baseline correction
    cfg             = [];
    cfg.baseline    = [-0.2 0];
    timelock{i}     = ft_timelockbaseline(cfg, timelock{i});
    
end 

% plot the ERPs over all sensors
figure 
for i = 1:blocktrials
    
    subplot(2,2,i);
    ft_singleplotER([], timelock{i});
    legend('averaged ERP over all sensors')
 
end

% save plots

% plot in interactive mode

% compute contrasts (would need to split the trials in conditions; i.e.
% easy vs difficult)


%% run statistics on ERPs





