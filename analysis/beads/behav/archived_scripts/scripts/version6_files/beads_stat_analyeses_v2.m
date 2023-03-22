% This script runs in combination with beads_prepro_v6.m 

% created: 20/12/2022
% modified: 21/12/2022 removed one sequence from sub 36 as the number of
% draws and the number of epochs for one sequence is inconsistent (biosemi
% did not record the events for 5 first draws in trial 14 (0.6 condition)).
% In the end decided to not include this trial in the EEG data and in
% model-fit data

% modified: 2/01/2023 added code to run ANOVAs with sites and laterality as
% factors and to multicompare (if needed).

% the script is used to:
% 1. Run rm anovas using the anovan.m function 
% 2. Compute differences in model_fit action values [draw again - max(urn choice)]
% 3. Extract cropped neural responses
% 3. Run regressions 
% 4. Run one sample t-tests 


%% LOAD ALL DATA %%

clear all 
clc

addpath(genpath('stored_workspaces'));

% all the data & outputs from preprocesing, io, and model-fiiting should be already saved in the working directory (see
% beads_prepro_v6.m)
% load ('workspace.mat')
load('two_param_fit.mat')
load('sub_data.mat')
load('draws_info')

%% DEFINE INIT PATHS & VARIABLES %%

% define vars
nsubs           = size(allsub_acc,1);
nconds          = size(allsub_acc,2);
totaltrials     = 52;
condtrials      = 26;

% define paths
wd              = pwd;
eeg_path        = fullfile(wd, 'cropped2');
scripts_path    = fullfile(wd, 'matlab_scripts');
addpath(genpath(eeg_path)); addpath(genpath(scripts_path))

%% RUN BEHAV STATISTICS - RM ANOVAs %%

% 1. 2x2 MIXED ANOVA on draws
% 2. 2x2 MIXED ANOVA on accuracy

% first add all the required data in one matrix 
subvec              = repmat(1:nsubs,1,4)';                             % create a vector with 4 copies participant number 
agentvec            = repmat([ones(1,nsubs*2) ones(1,nsubs*2)*2],1,1)'; % create a vector with 2 copies of agent type (indexed as 1=human, 2=io)
probvec             = repmat([ones(1,nsubs) ones(1,nsubs)*2],1,2)';     % create a vector with 2 copies of probability type (indexed as 1=0.8, 2=0.6)

% create 1 vec with all draws (human, io) 
drawsmat(:,1:2)     = allsub_draws;
drawsmat(:,3:4)     = allsubs_biodraws;
drawsvec            = drawsmat(:);

% create 1 vec with all acc (human, io) 
accmat(:,1:2)       = allsub_acc;
accmat(:,3:4)       = allsub_bioacc;
accvec              = accmat(:);

% run rm 2x2 anova on draws 
[pvals,~,stats] = anovan(drawsvec, {subvec agentvec probvec}, ... 
'model','interaction', 'random',1,'varnames',{'subvec' 'agentvec' 'probvec'})

% run rm 2x2 anova on accuracy 
[pvals,~,stats] = anovan(accvec, {subvec agentvec probvec}, ... 
'model','interaction', 'random',1,'varnames',{'subvec' 'agentvec' 'probvec'})

%% COMPUTE AQ DIFFERENCES & EXTRACT CROPPED NEURAL RESPONSES %%

for sub = 1:nsubs
    
    % extract model-fit Q values for this subject
    sub_modelQs     = modelfit_AQs{1,sub};
    sub_drawinfo    = allsub_drawinfo{1,sub};
    
    % first check model Q values and compute difference between max-urn
    % value and value for drawing again
    AQdiffs         = computeAQdiff_v2(sub_modelQs,condtrials,sub);
    
    % load subject-specific cropped MEEG file
    sub_eeg         = load(fullfile(eeg_path, sprintf('erpcropped_data_sub_%02d.mat', sub)));
    % sub_eeg         = load(fullfile(eeg_path, sprintf('tfrcropped_data_sub_%02d.mat', sub)));
    
    % run function to separate the "channels x samples x epochs" data based
    % on sites of interest (e.g., left and right frontal, parietal, etc..) and average by epochs to get one averaged
    % data point for each epoch/draw
    [allfrontal, allfc, allcp, allpar]  = eegExtract_v2(sub_eeg,sub_drawinfo);
    
    % we are done with this part!!
    
    %% RUN LINEAR REGRESSIONS AT THE SUBJECT LEVEL
    
    % run linear regressions and store the beta values for each laterality, site and
    % participant for ttests and ANOVAS
    % use the AQdiff values to predict neural activitity
    % [fbetas, fcbetas, cpbetas, pbetas] = runReg(AQdiffs, allfrontal, allfc, allcp, allpar);
    [fbetas, fcbetas, cpbetas, pbetas] = runRegV2(AQdiffs, allfrontal, allfc, allcp, allpar);
    % store beta values in arrays for all subs
    % loop over laterality 
    for l = 1:3
        
        allsubs_fb{l}(sub,1)    = fbetas{l};    % frontal (all,left,right)
        allsubs_fcb{l}(sub,1)   = fcbetas{l};   % frontocentral (all,left,right)
        allsubs_cpb{l}(sub,1)   = cpbetas{l};   % centroparietal (all,left,right)
        allsubs_pb{l}(sub,1)    = pbetas{l};    % parietal (all,left,right)
        
    end

end % end of subjects loop

%% Run ttests and save p values and t-stats

for l = 1:3
    
    % run frontal site 
    [h,p,ci,stats]      = ttest(allsubs_fb{l});
    frontal_pval(1,l)   = p;
    frontal_tstat(1,l)  = stats.tstat;
    frontal_sd(1,l)     = stats.sd;
    
    % run frontocentral site
    [h,p,ci,stats]      = ttest(allsubs_fcb{l});
    fc_pval(1,l)        = p;
    fc_tstat(1,l)       = stats.tstat;
    fc_sd(1,l)          = stats.sd;
    
    % run centroparietal site
    [h,p,ci,stats]      = ttest(allsubs_cpb{l});
    cp_pval(1,l)        = p;
    cp_tstat(1,l)       = stats.tstat;
    cp_sd(1,l)          = stats.sd;
    
    % run parietal site
    [h,p,ci,stats]      = ttest(allsubs_pb{l});
    p_pval(1,l)         = p;
    p_tstat(1,l)        = stats.tstat;
    p_sd(1,l)           = stats.sd;
    
end % end of laterality loop

%% ANOVAs to test sites and laterality (using the beta values)

% run anovas (mainly for multiple comparisons) using laterality (left,
% right) and site (frontal, fc, cp, par) as factors on the betas obtained
% at the subjects level

% first define factors 
% sites (frontal, fc, cp, par), and laterality (left,right)
sites           = repmat([ones(1,nsubs*2) ones(1,nsubs*2)*2 ones(1,nsubs*2)*3 ones(1,nsubs*2)*4],1,1)';
hemisph         = repmat([ones(1,nsubs) ones(1,nsubs)*2],1,1)';
laterality      = repmat(hemisph,4,1);

counter         = 0;
frontalr        = nan(40,1);
fcentralr       = nan(40,1);
cparr           = nan(40,1);
parr            = nan(40,1);

% extract arrays with the beta values from cells and concatinate them 
for l = 2:3
    
%     if l == 1
%         continue % we don't need the first as it includes left and right electrodes
%     end
    counter                 = counter + 1;
    frontalr(:,counter)     = allsubs_fb{l};
    fcentralr(:,counter)    = allsubs_fcb{l};
    cparr(:,counter)        = allsubs_cpb{l};
    parr(:,counter)         = allsubs_pb{l};
    
end % end of laterality cell

% concantinate the betas 
flater                      = frontalr(:);
fclater                     = fcentralr(:);
cplater                     = cparr(:);
plater                      = parr(:);

% add all sites in one array
allsites                    = cat(1, flater,fclater,cplater,plater);
subvec                      = repmat(1:nsubs,1,8)';                             % create a vector with 4 copies participant number 

% run rm anova 2 
stats                       = rm_anova2(allsites,subvec,laterality,sites,{'laterality' 'sites'});
