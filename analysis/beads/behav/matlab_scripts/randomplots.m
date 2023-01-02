% make more plots for Beads model fitting 

% load data
load('one_param_fit.mat')

% one and two parameter fit ll's 
easy_lls(:,1) = allsubs_ll_mb(:,1);
diff_lls(:,1) = allsubs_ll_mb(:,2);

% load 
load('two_param_fit.mat')

% one and two parameter fit ll's 
easy_lls(:,2) = allsubs_ll_mb(:,1);
diff_lls(:,2) = allsubs_ll_mb(:,2);

% plot participant average darws and model fit draws 
allsubjects = 1:40;

% all easy
easy_d(:,1) = modelfit_draws(:,1);
easy_d(:,2) = act_modelfit_draws(:,1);

% all difficult
diff_d(:,1) = modelfit_draws(:,2);
diff_d(:,2) = act_modelfit_draws(:,2);


%% Plot draws

% for easy and for difficult 
subplot(2,1,1)
plot(easy_d)
title('subjects vs model-fit draws 0.8')
xlabel('subjects')
ylabel('average draws')
legend('humans', 'model')

subplot(2,1,2)
plot(diff_d)
title('subjects vs model-fit draws 0.6')
xlabel('subjects')
ylabel('average draws')
legend('humans', 'model')

%% Plot ll

% for easy and for difficult 
subplot(2,1,1)
plot(easy_lls)
title('loglik for 0.8 trials')
xlabel('subjects')
ylabel('negative ll')
legend('one free param', 'two free params')

subplot(2,1,2)
plot(diff_lls)
title('loglik for 0.6 trials')
xlabel('subjects')
ylabel('negative ll')
legend('one free param', 'two free params')


%% PLOT MODEL FIT RESULTS %%

% allsub_draws(:,1)   = easy_avdraws;
% allsub_draws(:,2)   = diff_avdraws;
% allsub_acc(:,1)     = easy_avacc;
% allsub_acc(:,2)     = diff_avacc;
% % 
% % make plots
% makePlots(allsubs_beta_mb, allsubs_cs_mb, allsubs_ll_mb, allsub_draws)
% % makePlots(allsubs_cs_mb, allsubs_ll_mb, allsub_draws)
% 
% %% SAVE NEEDED FILES FOR ANALYSES %%
% 
% % save files that needed for statistical analyses
% save workspace
% save('sub_data', 'allsub_draws','allsub_acc', 'allsubs_biodraws', 'allsub_bioacc') % save subject and io outputs
% save('mf_data', 'modelfit_AQs', 'modelfit_draws', 'allsubs_cs_mb', 'allsubs_beta_mb', 'allsubs_ll_mb')
% 
% % lastly extract number of draws per subject 
% for sub =  1:nsubs
%     
%     subtemp         = allsub_alldata{1,sub};
%     subtotal(sub)   = sum(subtemp(:,5)+1); % +1 to account for last draw (which is the urn choice)
%     
% end
% 
% % transpose subtotal and save
% subtotal            = subtotal'; save('subtotal','subtotal', 'cond_data', 'allsub_alldata');

%% store model fit

save('two_param_fit.mat', 'act_modelfit_draws', 'modelfit_draws', 'actual_mdraws', 'mdraws', 'allsubs_cs_mb', 'allsubs_ll_mb', 'allsubs_AQs_mb', 'modelfit_AQs', 'allsubs_beta_mb')

