%% BEADS FORMAL PREPROCESSING AND ANALYSES SCRIPT W/ FIELDTRIP

% christinadelta 
% Date: Feb 2022
% last Update 25/02/2022

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
        
        % split the data into all-but-last samples (sample-again) and last
        % sample (price choice)
        
        
        
        
        
        
    end % end of blocks loop
    
    
end % end of subjects loop
