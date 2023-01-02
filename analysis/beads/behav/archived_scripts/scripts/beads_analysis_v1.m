
%% load workspacee 
load workspace 

% the next steps of analysis include dealing with the cropped EEG data, 
% running the association related model and running regression analysis.
% The association model will run similarly to the ideal observer, with the
% only difference that we will use cost to sample estimated in the
% parameterised model and not 0 or -0.25 as in the ideal observer. 

%% PREPARE FOR REGRESSION ANALYSIS OF CROPPED EEG WITH dAQ %%

for sub = 1:nsubs 
    
%     if sub == 9 | sub == 10 | sub == 23 | sub == 30 | sub == 34 | sub == 36 
%         continue
%     end

    if sub == 9 | sub == 10 | sub == 36 
        continue
    end

    
    % get the sum of the total draws (should be the same with the total number of
    % EEG epochs)
    sub_alldata             = allsub_alldata{1,sub};
    totaldraws              = sum(sub_alldata(:,5));
    
    % load the croped eeg data
    eeg_tmp                 = load(fullfile(croppedpath, sprintf('cropped_data_sub_%02d.mat',sub)));
    eeg_struct              = eeg_tmp.cropped_struct;
    eeg_data                = eeg_struct.data;
    eeg_conds               = eeg_struct.conds;
    eeg_events              = eeg_struct.events;
    
    % split the eeg data in frontal and parietal, fc, cp
    frontal_eeg             = eeg_data([1,2,3,4,17,18,19,20,21],:,:);
    fc_eeg                  = eeg_data([5,6,7,22,23,24,25],:,:);
    parietal_eeg            = eeg_data([11,12,13,14,15,29,30,31,32],:,:);
    cp_eeg                  = eeg_data([8,9,10,16,26,27,28],:,:);
    
    % combine frontal and fc and parietal and cp
    allfrontal_eeg = eeg_data([1,2,3,4,17,18,19,20,21,5,6,7,22,23,24,25],:,:);
    allparietal_eeg = eeg_data([11,12,13,14,15,29,30,31,32,8,9,10,16,26,27,28],:,:);
    
    
    % if length of total draws and events is not equal, correct it
    if totaldraws ~= length(eeg_events)
        
        totaldraws = length(eeg_events);
    end
    
   % average channels and samples for each trial and eeg matrix
    for trl = 1:totaldraws

        tmp_frontal                     = frontal_eeg(:,:,trl);
        tmp_parietal                    = parietal_eeg(:,:,trl);
        tmp_fc                          = fc_eeg(:,:,trl);
        tmp_cp                          = cp_eeg(:,:,trl);
        frontal_aveeg{1,sub}(trl,1)     = mean(tmp_frontal(:));  % average channels*samples for this trial
        parietal_aveeg{1,sub}(trl,1)    = mean(tmp_parietal(:)); % average channels*samples for this trial
        fc_aveeg{1,sub}(trl,1)          = mean(tmp_fc(:)); 
        cp_aveeg{1,sub}(trl,1)          = mean(tmp_cp(:));
        
        tmp_allfrontal                  = allfrontal_eeg(:,:,trl);
        tmp_allparietal                 = allparietal_eeg(:,:,trl);
        allfrontal_aveeg{1,sub}(trl,1)  = mean(tmp_allfrontal(:)); 
        allparietal_aveeg{1,sub}(trl,1) = mean(tmp_allparietal(:)); 
        
        

        clear tmp_frontal tmp_parietal tmp_fc tmp_cp tmp_allfrontal tmp_allparietal
    
    end % end of trl 
    
    % create a column of condition index to add in the eeg vector
    for i = 1:length(eeg_events)
        
        if eeg_events(i,1) == 1 | eeg_events(i,1) == 2
            eeg_events(i,2) = 1;
        elseif eeg_events(i,1) == 3 | eeg_events(i,1) == 4
            eeg_events(i,2) = 2;
        end
    end
    
    % create a column of trial/sequence number for each draw/epoch
    count = 0;
    subdraws = sub_alldata(:,5);
    
    for j = 1:totaltrials
        
        tmp_draws = subdraws(j,1);
        
        for d = 1:tmp_draws
            
            eeg_events(count+d,3) = j;

        end
        count = count + tmp_draws;
        
    end % end of j loop
    
     if length(eeg_events) > totaldraws
        
        diffr = length(eeg_events) - totaldraws;
        
        eeg_events(end,:) = [];
        
     end
    
    % short eeg points and events based on condition
    eeg_events(:,4)     = frontal_aveeg{1,sub};
    eeg_events(:,5)     = parietal_aveeg{1,sub};
    eeg_events(:,6)     = fc_aveeg{1,sub};
    eeg_events(:,7)     = cp_aveeg{1,sub};
    eeg_events(:,8)     = allfrontal_aveeg{1,sub};
    eeg_events(:,9)     = allparietal_aveeg{1,sub};
    eeg_sorted          = sortrows(eeg_events,2);

    %%% ok now the cropped eeg data are ready! 
    % run the association relate model 
    % this model will run similarly to the ideal observer. We will use the
    % Cs param estimated in model fitting but the rest will be as in the
    % ideal observer. Then the difference in AQ values will be calculated
    % as the difference between the max of the two urn choice values and
    % the draw again value. 

    % first add modelpath to the path
    addpath(genpath(iobserverpath));
    
    % extract this_sub estimated Cs from model fitting output
    sub_model = allsubs_model{1,sub};

    % define parameters of the model
    alpha       = 1;            % softmax stochasticity parameter (for fitting to human behaviour)
    Cw          = -10;        % cost for being wrong
    Cc          = 10;
    prob        = [0.8 0.6];    % proportion of the majority value in sequence (60:40 split in this case)
    cnt         = 1;            % counter
    count       = 1;            % second counter
    t           = 0;
    
    for cond = 1:conditions
        
        thiscond_data    = cond_data{1,subI}{1,cond};
        thiscond_seq     = allsequences{1,subI}{1,cond};

        % what is the probability of this cond? 
        if cond == 1
            thisq       = prob(1);
        else 
            thisq       = prob(2);
        end

        Cs         = sub_model(cond).params;
        
        % run association model
        [ll, pickTrial, dQvec, ddec, aQvec choice]  = estimateLikelihoodf_io(alpha,Cw,Cc,thisq,Cs,thiscond_seq,1);

        % extract AQ values for each condition
        cond_aQ                     = aQvec;
        
        % extract AQ values for each sequence 
        for trial = 1:length(cond_aQ)
            
            trial_aQ                = cond_aQ{1,trial};
            model_aQlen(count,1)    = length(trial_aQ);
            count                   = count+1; % update counter
            
            % loop over draws (AQs for that trial) 
            for draw = 1:size(trial_aQ,1)

                thisAQ                  = trial_aQ(draw,:);
                max_urn                 = max(thisAQ(1:2));
                allsub_dAQs(cnt,1)      = max_urn - thisAQ(3);
                allsub_dAQs(cnt,2)      = cond;
                allsub_dAQs(cnt,3)      = trial+t;
%                 allsub_drawAQs(cnt,1)   = trial_aQ(draw,3);
%                 allsub_drawAQs(cnt,2)   = cond;
%                 allsub_drawAQs(cnt,3)   = trial+t;

                cnt                            = cnt+1; % update counter

            end % end of draws (AQs) loop   
        end % end of trial loop
        
        t = t + (totaltrials/2); % update t
    end % end of conditions loop
    
    % add new trial number to the sorted EEG matrix (that corresponds to the AQs trial number)
    for cond = 1:conditions
        
        thiscond_draws(:,cond) = cond_data{1,sub}{1,cond}(:,5);
        thiscond_draws = thiscond_draws(:);
   
    end % end of condition loop
    
    % add a new index for each draw     
     count = 0;
     for ii = 1:totaltrials
         
         tmp_draws = thiscond_draws(ii,1);
         for jj = 1:tmp_draws
             eeg_sorted(count+jj,10) = ii;
             
         end % end of jj loop
         count = count + tmp_draws;
     end % end of ii loop
    
    % the EEG vector and the AQs vector need to be of equal sizes. However
    % on every sequence/trial the model draws more than then participant.
    % Thus, for every sequence if the model AQs are more than EEG epochs,
    % fill the EEG vector with nans (and vice versa).
    
    count = 1;
    cnt = 1;
    for trl = 1:totaltrials 
        
        tmp_f           = find(eeg_sorted(:,10)==trl);
        tmp_frontal     = eeg_sorted((tmp_f),4);
        tmp_parietal    = eeg_sorted((tmp_f),5);
        tmp_fc          = eeg_sorted((tmp_f),6);
        tmp_cp          = eeg_sorted((tmp_f),7);
        tmp_allfrontal  = eeg_sorted((tmp_f),8);
        tmp_allparietal = eeg_sorted((tmp_f),9);
        
        
        tmp_aq      = find(allsub_dAQs(:,3)==trl);
        tmp_daq     = allsub_dAQs((tmp_aq),1);
%         tmp_diffaq  = allsub_diffAQs((tmp_aq),1);
%         tmp_drawaq  = allsub_drawAQs((tmp_aq),1);
        
        if length(tmp_aq) > length(tmp_f)
            
            s = length(tmp_f)+1;
            e = length(tmp_aq);
            tmp_frontal(s:e,1) = nan;
            tmp_parietal(s:e,1) = nan;
            tmp_fc(s:e,1) = nan;
            tmp_cp(s:e,1) = nan;
            tmp_allfrontal(s:e,1) = nan;
            tmp_allparietal(s:e,1) = nan;

                
        elseif length(tmp_aq) < length(tmp_f)
            
            s = length(tmp_aq)+1;
            e = length(tmp_f);
            tmp_daq(s:e,1) = nan;
%             tmp_diffaq(s:e,1) = nan;
%             tmp_drawaq(s:e,1) = nan;

        end
        
        % add them all in one vector 
        for jj = 1:length(tmp_frontal)

            eeg_f(count,1) = tmp_frontal(jj,1);
            eeg_p(count,1) = tmp_parietal(jj,1);
            eeg_fc(count,1) = tmp_fc(jj,1);
            eeg_cp(count,1) = tmp_cp(jj,1);
            eeg_allf(count,1) = tmp_allfrontal(jj,1);
            eeg_allp(count,1) = tmp_allparietal(jj,1);

            count = count+1;

        end
        
        for ii = 1:length(tmp_daq)
                
            dAQs(cnt,1) = tmp_daq(ii,1);
            cnt = cnt+1;

        end
        
        clear tmp_frontal tmp_parietal tmp_daq tmp_diffaq tmp_drawaq tmp_f tmp_aq tmp_fc tmp_cp tmp_allfrontal 
        clear tmp_allparietal

    end % end of trial loop
    
    %% RUN REGRESSIONS
    
    % Parietal Regressions
    
    % regression 1 - dAQs, eeg_p
    mdl = fitlm(dAQs,eeg_p)
    p_betas_dAQS(sub,1) = table2array(mdl.Coefficients(2,1));
    clear mdl
    
    % regression 2 - dAQs, eeg_cp
    mdl = fitlm(dAQs,eeg_cp)
    cp_betas_dAQS(sub,1) = table2array(mdl.Coefficients(2,1));
    clear mdl
    
    % regression 3 - dAQs, eeg_allp
    mdl = fitlm(dAQs,eeg_allp)
    allp_betas_dAQS(sub,1) = table2array(mdl.Coefficients(2,1));
    clear mdl
    
    % Frontal Regressions 
    
    % regression 4 - dAQs, eeg_f
    mdl = fitlm(dAQs,eeg_f)
    f_betas_dAQS(sub,1) = table2array(mdl.Coefficients(2,1));
    clear mdl
    
    % regression 5 - dAQs, eeg_fc
    mdl = fitlm(dAQs,eeg_fc)
    fc_betas_dAQS(sub,1) = table2array(mdl.Coefficients(2,1));
    clear mdl
    
    % regression 6 - dAQs, eeg_allf
    mdl = fitlm(dAQs,eeg_allf)
    allf_betas_dAQS(sub,1) = table2array(mdl.Coefficients(2,1));
    clear mdl
    
    
    %% RUN PEARSON'S CORRELATIONS AT SUBJECT LEVEL
    
    % for each subject run pearson's r 
    % parietal correlations
    % parietal & diff in AQs
    R = corrcoef(dAQs, eeg_p, 'Rows','complete')
    % z = atanh(R); % z-transform
    all_parietal_r(sub,1) = R(1,2);
    
    clear R 
    
    % centroparietal & diff in AQs
    R = corrcoef(dAQs, eeg_cp, 'Rows','complete')
    % z = atanh(R); % z-transform
    all_cp_r(sub,1) = R(1,2);
     clear R 
     
    % parietal + centroparietal & diff in AQs
    R = corrcoef(dAQs, eeg_allp, 'Rows','complete')
    % z = atanh(R); % z-transform
    all_cp_par_r(sub,1) = R(1,2);
    clear R 
    
    % frontal correlations
    % frontal & diff in AQs
    R = corrcoef(dAQs, eeg_f, 'Rows','complete')
    % z = atanh(R); % z-transform
    all_frontal_r(sub,1) = R(1,2);
    
    clear R 
    
    % fc & diff in AQs
    R = corrcoef(dAQs, eeg_fc, 'Rows','complete')
    % z = atanh(R); % z-transform
    all_fc_r(sub,1) = R(1,2);
     clear R 
     
    % frontal + fc & diff in AQs
    R = corrcoef(dAQs, eeg_allf, 'Rows','complete')
    % z = atanh(R); % z-transform
    all_fc_front_r(sub,1) = R(1,2);
    clear R 
    
    clear sub_alldata totaldraws eeg_tmp eeg_struct eeg_data eeg_conds eeg_events frontal_eeg parietal_eeg dQvec
    
    clear tmp_frontal tmp_parietal sub_draws count tmp_draws eeg_sorted thiscond_data thiscond_seq Cs thisq ll pickTrial 
    clear ddec aQvec choice cond_aQ trial_aQ model_aQlen count thisAQ max_urn allsub_dAQs allsub_diffAQs allsub_drawAQs t
    
    clear thiscond_draws dAQs eeg_f eeg_fc eeg_allf eeg_p eeg_cp eeg_allp
    
    
end % end of subjects loop

%% Run one-sample t-tests

% Run 2nd level stats on regressions
% parietal ttests
[h,p,ci,stats] = ttest(p_betas_dAQS)
[h,p,ci,stats] = ttest(cp_betas_dAQS)
[h,p,ci,stats] = ttest(allp_betas_dAQS)

% frontal ttests
[h,p,ci,stats] = ttest(f_betas_dAQS)
[h,p,ci,stats] = ttest(fc_betas_dAQS)
[h,p,ci,stats] = ttest(allf_betas_dAQS)

%%%%%%%
%%%%%%%

% z-transform corr coefficients
all_parietal_z = atanh(all_parietal_r); % z-transform
all_cp_z = atanh(all_cp_r); %
all_cp_par_z = atanh(all_cp_par_r);

all_frontal_z = atanh(all_frontal_r);
all_fc_z = atanh(all_fc_r);
all_fc_front_z = atanh(all_fc_front_r);

% Run 2nd level stats on correlations 
% parietal sites
[h,p,ci,stats] = ttest(all_parietal_z)
[h,p,ci,stats] = ttest(all_cp_z)
[h,p,ci,stats] = ttest(all_cp_par_z)

% frontal ttests
[h,p,ci,stats] = ttest(all_frontal_z)
[h,p,ci,stats] = ttest(all_fc_z)
[h,p,ci,stats] = ttest(all_fc_front_z)
