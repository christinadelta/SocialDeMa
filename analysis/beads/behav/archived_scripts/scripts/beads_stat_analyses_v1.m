% This script runs in combination with beads_prepro_v6.m 

% created: 14/11/22

% modified: 16:12/2022
% This version uses the model AQ values up until the moment that the model
% should have stopped (not the values up until participants stopping
% point), thus, some of the EEG epochs are removed to bring the model AQ
% values and participants' EEG epochs at the same length 

% Version 2 uses all model Q values (up until participants stopping
% points) and all EEG epochs (without removing anything)


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
% subvec              = repmat(1:nsubs,1,4)';                             % create a vector with 4 copies participant number 
% agentvec            = repmat([ones(1,nsubs*2) ones(1,nsubs*2)*2],1,1)'; % create a vector with 2 copies of agent type (indexed as 1=human, 2=io)
% probvec             = repmat([ones(1,nsubs) ones(1,nsubs)*2],1,2)';     % create a vector with 2 copies of probability type (indexed as 1=0.8, 2=0.6)
% 
% % create 1 vec with all draws (human, io) 
% drawsmat(:,1:2)     = allsub_draws;
% drawsmat(:,3:4)     = allsubs_biodraws;
% drawsvec            = drawsmat(:);
% 
% % create 1 vec with all acc (human, io) 
% accmat(:,1:2)       = allsub_acc;
% accmat(:,3:4)       = allsub_bioacc;
% accvec              = accmat(:);
% 
% % run rm 2x2 anova on draws 
% [pvals,~,stats] = anovan(drawsvec, {subvec agentvec probvec}, ... 
% 'model','interaction', 'random',1,'varnames',{'subvec' 'agentvec' 'probvec'})
% 
% % run rm 2x2 anova on accuracy 
% [pvals,~,stats] = anovan(accvec, {subvec agentvec probvec}, ... 
% 'model','interaction', 'random',1,'varnames',{'subvec' 'agentvec' 'probvec'})

%% COMPUTE AQ DIFFERENCES & EXTRACT CROPPED NEURAL RESPONSES %%

for sub = 1:nsubs
    
%     if sub == 36 % don't include sub 36 for now
%         continue 
%     end
%     
    % extract model-fit Q values for this subject
    sub_modelQs     = modelfit_AQs{1,sub};
    sub_drawinfo    = allsub_drawinfo{1,sub};
    
    % load subject-specific cropped MEEG file
    sub_eeg     = load(fullfile(eeg_path, sprintf('erpcropped_data_sub_%02d.mat', sub)));
    % sub_eeg         = load(fullfile(eeg_path, sprintf('tfrcropped_data_sub_%02d.mat', sub)));
    % subdraw         = subtotal(sub);
    
    % run function to separate the "channels x samples x epochs" data based
    % on channels of interest and average by epochs to get one averaged
    % data point for each epoch/draw
    [allfrontal, allfc, allcp, allpar]  = eegExtract(sub_eeg, sub_drawinfo, sub);
    
    for cond = 1:nconds
        
        % run function to compute differences in Q values
        cond_modelQs                    = sub_modelQs{1,cond};
        AQdiffs                         = computeAQdiff(cond_modelQs, condtrials);
        all_AQdiffs{1,cond}             = AQdiffs; % just save the AQ difference values; will need them later

    end
    
    % make sure that AQdiff values and eeg data-points are of the same
    % length for each trial
    % [new_AQs, f, fc, cp, par]           = checkValues(allfrontal, allfc, allcp, allpar, all_AQdiffs);
    [new_AQs, f, fc, cp, par]           = checkValuesV2(allfrontal, allfc, allcp, allpar, all_AQdiffs, sub);
    
    %% RUN LINEAR REGRESSIONS AT THE SUBJECT LEVEL
    
    % use the AQdiff values to predict neural activitity
    % [fbetas, fcbetas, cpbetas, pbetas] = runReg(new_AQs, f, fc, cp, par);
    [fbetas, fcbetas, cpbetas, pbetas] = runRegV2(new_AQs, f, fc, cp, par);
    
    % add participant beta values (for each regression in one array
    % frontal
    allsub_fbetas(sub,1) = fbetas.f;
    allsub_lfbetas(sub,1) = fbetas.lf;
    allsub_rfbetas(sub,1) = fbetas.rf;
    
    % frontocentral
    allsub_fcbetas(sub,1) = fcbetas.fc;
    allsub_lfcbetas(sub,1) = fcbetas.lfc;
    allsub_rfcbetas(sub,1) = fcbetas.rfc;
    
    % centroparietal
    allsub_cpbetas(sub,1) = cpbetas.cp;
    allsub_lcpbetas(sub,1) = cpbetas.lcp;
    allsub_rcpbetas(sub,1) = cpbetas.rcp;
    
    % parietal 
    allsub_pbetas(sub,1) = pbetas.p;
    allsub_lpbetas(sub,1) = pbetas.lp;
    allsub_rpbetas(sub,1) = pbetas.rp;
    
end % end of subject loop

%% Run one-sample ttests (v1)

% % frontal sites ttests
% [h,p,ci,stats] = ttest(allsub_fbetas)
% [h,p,ci,stats] = ttest(allsub_lfbetas)
% [h,p,ci,stats] = ttest(allsub_rfbetas)
% 
% % frontocentral sites ttests
% [h,p,ci,stats] = ttest(allsub_fcbetas)
% [h,p,ci,stats] = ttest(allsub_lfcbetas)
% [h,p,ci,stats] = ttest(allsub_rfcbetas)
% 
% % centroparietal sites ttests
% [h,p,ci,stats] = ttest(allsub_cpbetas)
% [h,p,ci,stats] = ttest(allsub_lcpbetas)
% [h,p,ci,stats] = ttest(allsub_rcpbetas)
% 
% % parietal sites ttests
% [h,p,ci,stats] = ttest(allsub_pbetas)
% [h,p,ci,stats] = ttest(allsub_lpbetas)
% [h,p,ci,stats] = ttest(allsub_rpbetas)


%% Run one-sample ttests (v2)

% convert tables to arrays 
% frontal
all_fbetas = table2array(allsub_fbetas);
all_lfbetas = table2array(allsub_lfbetas);
all_rfbetas = table2array(allsub_rfbetas);

% fc
all_fcbetas = table2array(allsub_fcbetas);
all_lfcbetas = table2array(allsub_lfcbetas);
all_rfcbetas = table2array(allsub_rfcbetas);

% cp
all_cpbetas = table2array(allsub_cpbetas);
all_lcpbetas = table2array(allsub_lcpbetas);
all_rcpbetas = table2array(allsub_rcpbetas);

% par
all_pbetas = table2array(allsub_pbetas);
all_lpbetas = table2array(allsub_lpbetas);
all_rpbetas = table2array(allsub_rpbetas);

% this analysis is at group level 
% run one-sample ttests using the betas (from first level) to test that the samples differ
% from zero

% frontal sites ttests
[h,p,ci,stats] = ttest(all_fbetas)
[h,p,ci,stats] = ttest(all_lfbetas)
[h,p,ci,stats] = ttest(all_rfbetas)

% frontocentral sites ttests
[h,p,ci,stats] = ttest(all_fcbetas)
[h,p,ci,stats] = ttest(all_lfcbetas)
[h,p,ci,stats] = ttest(all_rfcbetas)

% centroparietal sites ttests
[h,p,ci,stats] = ttest(all_cpbetas)
[h,p,ci,stats] = ttest(all_lcpbetas)
[h,p,ci,stats] = ttest(all_rcpbetas)

% parietal sites ttests
[h,p,ci,stats] = ttest(all_pbetas)
[h,p,ci,stats] = ttest(all_lpbetas)
[h,p,ci,stats] = ttest(all_rpbetas)