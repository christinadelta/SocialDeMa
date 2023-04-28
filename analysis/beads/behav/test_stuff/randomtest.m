%% MODEL FIT 1 %% 

% first add model fitting to the path
addpath(genpath(bmodelfitpath));


% RUN MODEL 1:
% FREE PARAMS: ONLY BETA 
% define a few parameters 
R.error             = -10;          % cost for being wrong
R.diff              = -20;          % The difference between the rewards for being correct (in this case no reward 10) and the cost of being wrong (-10).
R.correct           = 10;           % reward for being correct
R.q                 = [0.8 0.6];    % proportion of the majority value in sequence (60:40 split in this case)
R.initsample        = -0.25;        % the cost to sample
R.freeparams        = 1;            % only beta is used
range_betas         = 1;
R.modelnum          = 1;

for beta_model = 1:range_betas

    % init beta value
    % this_beta               = exprnd(10);
    this_beta               = 3;
    R.initbeta              = this_beta;
    beta_vals(beta_model)   = this_beta;

    for sub = 1:nsubs

        subdata         = cond_data{1,sub};         % all data matrix
        sub_choicesvec  = allsub_choiceVec{1,sub};  % subject choices (two 26 by 3 matricies)
        sub_sequence    = allsub_sequences{1,sub};  % sequence that this-subject was presented with

        for cond = 1:conditions

            % extract condition data
            sub_cond        = subdata{1,cond};
            cond_choices    = sub_choicesvec{1,cond};
            cond_sequence   = sub_sequence{1,cond};

            % what is the probability of this cond? 
            if cond == 1

                thisq       = R.q(1);
            else 
                thisq       = R.q(2);
            end
        
            R.thisq         = thisq;

            % extract sequences from cell and store in matrix (26x10)
            for s = 1:size(cond_sequence,2)
                thiscond_seqmat(s,:)                        = cond_sequence{1,s};
            end
            
            % what is the urntype?
            urntype                                         = sub_cond(:,4);

            for u = 1:length(urntype)
                if urntype(u) == 0 % if green urn switch index coding
                    seq_ones                                = find(thiscond_seqmat(u,:) == 1);
                    seq_twos                                = find(thiscond_seqmat(u,:) == 2);
                    thiscond_seqmat(u,seq_ones)             = 2;
                    thiscond_seqmat(u,seq_twos)             = 1;
                end 
            end

            % recode 2s to 0s for backward induction 
            thiscond_seqmat(find(thiscond_seqmat==2))       = 0;
        
            % fit free parameter using Bruno's version of the model
            [minParams, ll, Qsad, cprob, model_samples,model_urnchoice]     = bayesbeads_b(thiscond_seqmat, cond_choices, R);

            
            model1(cond).param                              = minParams;
            model1(cond).ll                                 = ll;
            model1(cond).Qsad                               = Qsad;
            model1(cond).cprob                              = cprob;
            model1(cond).samples                            = model_samples;

            % all_model_samples(sub,cond,beta_model)          = mean(model_samples);
            model1_samples(sub,cond)                        = mean(model_samples); % for ploting only right now
            allsubs_beta_model1(sub,cond)                 = minParams;
            allsubs_ll_model1(sub,cond)                     = ll;

            % compute model performance (to be used for statistical
            % analysis 
            for t = 1:size(model_urnchoice,1)
                
                if (model_urnchoice(t) == 1 & urntype(t) == 1) | (model_urnchoice(t) == 2 & urntype(t) == 0)
                    model_choice(t) = 1;
                else
                    model_choice(t) = 0;
                end % end of condition
            end % end of trials loop

            model1_acc(sub,cond)                            = mean(model_choice);
            model1(cond).acc                                = mean(model_choice);

        end % end of conditions loop

        % store subject and model data
        all_model1{beta_model,sub}                          = model1;

    end % end of subjects loop

    model_sampling_rate(beta_model,:)                       = mean(model1_samples,1); % average across subjects 

end % end of models loop

% choose the first one for ploting and comparing with model 2
% model1_1 = all_model_samples(:,:,1)
% bar(model_sampling_rate) % check sampling rates

clear minParams ll Qsad cprob model_samples model_urnchoice R










%% MODEL FIT 2 %%

% RUN MODEL 2:
% free parameters: [cost-to-sample]
% define a few parameters 
R.error             = -10;          % cost for being wrong
R.diff              = -20;          % The difference between the rewards for being correct (in this case no reward 10) and the cost of being wrong (-10).
R.correct           = 10;           % reward for being correct
R.q                 = [0.8 0.6];    % proportion of the majority value in sequence (60:40 split in this case)
R.initsample        = -0.25;        % the cost to sample
R.freeparams        = 1;            %  Cs and beta are used
R.initbeta          = 3;            % 
R.modelnum          = 2;

for sub = 1:nsubs

    subdata         = cond_data{1,sub};         % all data matrix
    sub_choicesvec  = allsub_choiceVec{1,sub};  % subject choices (two 26 by 3 matricies)
    sub_sequence    = allsub_sequences{1,sub};  % sequence that this-subject was presented with

    for cond = 1:conditions

        % extract condition data
        sub_cond        = subdata{1,cond};
        cond_choices    = sub_choicesvec{1,cond};
        cond_sequence   = sub_sequence{1,cond};

        % what is the probability of this cond? 
        if cond == 1

            thisq       = R.q(1);
        else 
            thisq       = R.q(2);
        end
    
        R.thisq         = thisq;

        % extract sequences from cell and store in matrix (26x10)
        for s = 1:size(cond_sequence,2)
            thiscond_seqmat(s,:)                        = cond_sequence{1,s};
        end
        
        % what is the urntype?
        urntype                                         = sub_cond(:,4);

        for u = 1:length(urntype)
            if urntype(u) == 0 % if green urn switch index coding
                seq_ones                                = find(thiscond_seqmat(u,:) == 1);
                seq_twos                                = find(thiscond_seqmat(u,:) == 2);
                thiscond_seqmat(u,seq_ones)             = 2;
                thiscond_seqmat(u,seq_twos)             = 1;
            end 
        end

        % recode 2s to 0s for backward induction 
        thiscond_seqmat(find(thiscond_seqmat==2))       = 0;
    
        % fit free parameter using Bruno's version of the model
        [minParams, ll, Qsad, cprob, model_samples,model_urnchoice]     = bayesbeads_b(thiscond_seqmat, cond_choices, R);

        
        model2(cond).param                              = minParams;
        model2(cond).ll                                 = ll;
        model2(cond).Qsad                               = Qsad;
        model2(cond).cprob                              = cprob;
        model2(cond).samples                            = model_samples;

        % all_model_samples(sub,cond,beta_model)          = mean(model_samples);
        model2_samples(sub,cond)                        = mean(model_samples); % for ploting only right now
        allsubs_Cs_model2(sub,cond)                     = minParams;
        allsubs_ll_model2(sub,cond)                     = ll;

        % compute model performance (to be used for statistical
        % analysis 
        for t = 1:size(model_urnchoice,1)
            
            if (model_urnchoice(t) == 1 & urntype(t) == 1) | (model_urnchoice(t) == 2 & urntype(t) == 0)
                model_choice(t) = 1;
            else
                model_choice(t) = 0;
            end % end of condition
        end % end of trials loop

        % compute accuracy
        model2_acc(sub,cond)                            = mean(model_choice);
        model2(cond).acc                                = mean(model_choice);

    end % end of conditions loop

    % store subject and model data
    all_model2{1,sub}                                   = model2;

end % end of subjects loop

clear minParams ll Qsad cprob model_samples model_urnchoice R

%% RUN ANOVAS AND MULRICOMPARE %%

% add draws and performance in one vec
all_acc             = [easy_avacc diff_avacc];
all_draws           = [easy_avdraws diff_avdraws];

% make a struct with all the vectors needed for the analysis
anova_struct        = struct('all_draws', all_draws, 'all_acc', all_acc,...
    'all_ioacc', all_ioacc, 'all_iodraws', all_iodraws, 'model2_samples', model2_samples, 'model2_acc', model2_acc);

% run anovas and pairwise comparisons 
output_struct_one   = runBehavStats(nsubs, anova_struct); % output will be used for plotting 

%% MODEL FIT 3 %%

% RUN MODEL 2:
% free parameters: [beta, cost-to-sample]
% define a few parameters 
R.error             = -10;          % cost for being wrong
R.diff              = -20;          % The difference between the rewards for being correct (in this case no reward 10) and the cost of being wrong (-10).
R.correct           = 10;           % reward for being correct
R.q                 = [0.8 0.6];    % proportion of the majority value in sequence (60:40 split in this case)
R.initsample        = -0.25;        % the cost to sample
R.freeparams        = 2;            %  Cs and beta are used
R.initbeta          = 3;            % 


for sub = 1:nsubs

    subdata         = cond_data{1,sub};         % all data matrix
    sub_choicesvec  = allsub_choiceVec{1,sub};  % subject choices (two 26 by 3 matricies)
    sub_sequence    = allsub_sequences{1,sub};  % sequence that this-subject was presented with

    for cond = 1:conditions

        % extract condition data
        sub_cond        = subdata{1,cond};
        cond_choices    = sub_choicesvec{1,cond};
        cond_sequence   = sub_sequence{1,cond};

        % what is the probability of this cond? 
        if cond == 1

            thisq       = R.q(1);
        else 
            thisq       = R.q(2);
        end
    
        R.thisq         = thisq;

        % extract sequences from cell and store in matrix (26x10)
        for s = 1:size(cond_sequence,2)
            thiscond_seqmat(s,:)                        = cond_sequence{1,s};
        end
        
        % what is the urntype?
        urntype                                         = sub_cond(:,4);

        for u = 1:length(urntype)
            if urntype(u) == 0 % if green urn switch index coding
                seq_ones                                = find(thiscond_seqmat(u,:) == 1);
                seq_twos                                = find(thiscond_seqmat(u,:) == 2);
                thiscond_seqmat(u,seq_ones)             = 2;
                thiscond_seqmat(u,seq_twos)             = 1;
            end 
        end

        % recode 2s to 0s for backward induction 
        thiscond_seqmat(find(thiscond_seqmat==2))       = 0;
    
        % fit free parameter using Bruno's version of the model
        [minParams, ll, Qsad, cprob, model_samples,model_urnchoice]     = bayesbeads_b(thiscond_seqmat, cond_choices, R);

        
        model3(cond).param                              = minParams;
        model3(cond).ll                                 = ll;
        model3(cond).Qsad                               = Qsad;
        model3(cond).cprob                              = cprob;
        model3(cond).samples                            = model_samples;
        model3_samples(sub,cond)                        = mean(model_samples); % for ploting only right now
        allsubs_cs_model3(sub,cond)                     = minParams(1);
        allsubs_beta_model3(sub,cond)                   = minParams(2);
        allsubs_ll_model3(sub,cond)                     = ll;


        % compute model performance
        for t = 1:size(model_urnchoice,1)
                
            if (model_urnchoice(t) == 1 & urntype(t) == 1) | (model_urnchoice(t) == 2 & urntype(t) == 0)
                model_choice(t) = 1;
            else
                model_choice(t) = 0;
            end % end of condition
        end % end of trials loop

        model3_acc(sub,cond)                            = mean(model_choice);
        model3(cond).acc                                = mean(model_choice);

    end % end of conditions loop

    all_model3{1,sub}                                   = model2;

end % end of subjects loop

% model2_samples_mean = mean(model2_samples_plot,1)
% allmodels_barplot = [model_sampling_rate(1,:); model2_samples_mean]
% bar(allmodels_barplot)

clear minParams ll Qsad cprob model_samples model_urnchoice R

%% RUN ANOVAS & PAIRWISE COMPARISONS %%

% factors: Agent type (human, io, beta, beta_cs), probability 
anova_struct        = struct('all_draws', all_draws, 'all_acc', all_acc,'all_ioacc', all_ioacc,...
    'all_iodraws',all_iodraws,'model1_samples',model1_samples,'model1_acc',model1_acc,'model2_samples',...
    model2_samples,'model2_acc',model2_acc, 'model3_samples', model3_samples, 'model3_acc',model3_acc);

output_struct_two   = runBehavStats(nsubs, anova_struct); % output will be used for plotting 

%% PARAMETER RECOVERY FOR MODEL 1 (BETA) %%

% in parameter recover simulate 52 datasets and fit the model N
% times with a range of Beta parameter values (in model 1) and a range of
% Beta and Cs parameters in model 2. To see whether parameter recovery is 
% good I plot the simulated vs fitted parameter values in scatter plots.

% define variables for data simulation
simvars.ntrials     = 52;
simvars.maxDraws    = 10;
simvars.qvals       = [0.8 0.6];
simvars.conditions  = 2;

% define parameters for simulations
% simR.rangeCs        = [-5:0.25:0];
simR.correct            = 10;
simR.error              = -10;
simR.diff               = -20;
simR.initsample         = -0.25; % cost-to-sample is fixed in beta model
simR.freeparams         = 1; 
simR.modelnum           = 1;

% deal with betas
beta_bounds             =[10 15];        % maximum and minimum value of beta?
nbins                   = 10;
allbetas                = linspace(beta_bounds(1), beta_bounds(2), nbins+1);

nReps                   = length(allbetas);

% RUN MODEL 1 -- BETA
for rep = 1:nReps

    simR.initbeta           = allbetas(1,rep); % this rep beta value?
    
    % simulate dataset of 52 sequences/trials (26 easy and 26 difficult ones)
    simoutput = simBeadsData(simvars,simR);

    % store the main simulated parameters and samples
    simulated_model1(rep).samples          = simoutput.simdraws;
    simulated_model1(rep).avsamples        = simoutput.avsamples;
    
    for cond = 1:conditions

        cond_simsequence                                        = simoutput.simsequences{1,cond};
        cond_simvecs                                            = simoutput.simchoicevec{1,cond};
        simR.thisq                                              = simvars.qvals(cond);

        % what is the simulated urntype?
        simurntype                                              = simoutput.simurns(:,cond);

        for u = 1:length(simurntype)
            if simurntype(u) == 0 % if green urn switch index coding
                seq_ones                                    = find(cond_simsequence(u,:) == 1);
                seq_twos                                    = find(cond_simsequence(u,:) == 2);
                cond_simsequence(u,seq_ones)                = 2;
                cond_simsequence(u,seq_twos)                = 1;
            end 
        end

        % recode 2s to 0s for backward induction 
        cond_simsequence(find(cond_simsequence==2))         = 0;

        
        % now that we have simulated sequences and responses, fit the model
        % with the simulated R struct
        [minParams,ll,Qsad,cprob,fitted_samples,fitted_urnchoice] = bayesbeads_b(cond_simsequence, cond_simvecs, simR);

        % compute fitted/recovered performance (to be used for statistical
        % analysis 
        for t = 1:size(fitted_urnchoice,1)
            
            if (fitted_urnchoice(t) == 1 & simurntype(t,1) == 1) | (fitted_urnchoice(t) == 2 & simurntype(t,1) == 0)
                fitted_choice(t) = 1;
            else
                fitted_choice(t) = 0;
            end % end of condition
        end % end of trials loop
        
        performance(cond)                       = mean(fitted_choice);
        fitted_avsamples(cond)                  = mean(fitted_samples);
        
        % store the main simulated parameters and samples
        simulated_model1(rep).beta(cond)       = simR.initbeta;
        simulated_model1(rep).correct(cond)    = simoutput.performance(cond);

        % store fitted parameters and samples
        fitted_model1(rep).samples(:,cond)     = fitted_samples;
        fitted_model1(rep).beta(cond)          = minParams(1);
        fitted_model1(rep).correct(cond)       = performance(cond);
        fitted_model1(rep).avsamples(cond)     = fitted_avsamples(cond);
        fitted_model1(rep).ll(cond)            = ll;
        
    end 

end % end of model 1 repetitions loop

clear simR fitted_samples fitted_urnchoice fitted_choice

%% PARAMETER RECOVER FOR MODEL 2 (COST-SAMPLE) %%




%% PARAMETER RECOVER FOR MODEL 2 (BETA & COST-SAMPLE) %%

% define parameters for simulations
% simR.rangeCs        = [-5:0.25:0];
simR.correct            = 10;
simR.error              = -10;
simR.diff               = -20;
simR.freeparams         = 2; 

% deal with betas
beta_bounds             = [0.1 5];        % maximum and minimum value of beta?
cs_bounds               = [-2 0]; % maximum and minimum cost to sample
nbins                   = [10 8];
allbetas                = linspace(beta_bounds(1), beta_bounds(2), nbins(1)+1);
allCs                   = linspace(cs_bounds(1), cs_bounds(2), nbins(2)+1);

% RUN PARAMETER RECOVERY FORM MODEL 2
for cs = 1:length(allCs)

    simR.initsample = allCs(cs); % extract cost-sample for this loop

    for beta = 1:length(allbetas)

        simR.initbeta = allbetas(beta);

        % simulate dataset of 52 sequences/trials (26 easy and 26 difficult ones)
        simoutput = simBeadsData(simvars,simR);

        % store the main simulated parameters and samples
        simulated_model2(beta).samples          = simoutput.simdraws;
        simulated_model2(beta).avsamples        = simoutput.avsamples;

        for cond = 1:conditions

            cond_simsequence                                        = simoutput.simsequences{1,cond};
            cond_simvecs                                            = simoutput.simchoicevec{1,cond};
            simR.thisq                                              = simvars.qvals(cond);
    
            % what is the simulated urntype?
            simurntype                                         = simoutput.simurns(:,cond);
    
            for u = 1:length(simurntype)
                if simurntype(u) == 0 % if green urn switch index coding
                    seq_ones                                = find(cond_simsequence(u,:) == 1);
                    seq_twos                                = find(cond_simsequence(u,:) == 2);
                    cond_simsequence(u,seq_ones)             = 2;
                    cond_simsequence(u,seq_twos)             = 1;
                end 
            end
    
            % recode 2s to 0s for backward induction 
            cond_simsequence(find(cond_simsequence==2))       = 0;
    
            
            % now that we have simulated sequences and responses, fit the model
            % with the simulated R struct
            [minParams,ll,Qsad,cprob,fitted_samples,fitted_urnchoice] = bayesbeads_b(cond_simsequence, cond_simvecs, simR);
    
            % compute fitted/recovered performance (to be used for statistical
            % analysis 
            for t = 1:size(fitted_urnchoice,1)
                
                if (fitted_urnchoice(t) == 1 & simurntype(t,1) == 1) | (fitted_urnchoice(t) == 2 & simurntype(t,1) == 0)
                    fitted_choice(t) = 1;
                else
                    fitted_choice(t) = 0;
                end % end of condition
            end % end of trials loop

            performance(cond)                       = mean(fitted_choice);
            fitted_avsamples(cond)                  = mean(fitted_samples);
            
            % store the main simulated parameters and samples
            simulated_model2(beta).beta(cond)       = simR.initbeta;
            simulated_model2(beta).correct(cond)    = simoutput.performance(cond);
            simulated_model2(beta).cs(cond)         = simR.initsample;
    
            % store fitted parameters and samples
            fitted_model2(beta).samples(:,cond)     = fitted_samples;
            fitted_model2(beta).beta(cond)          = minParams(2);
            fitted_model2(beta).cs(cond)            = minParams(1);
            fitted_model2(beta).correct(cond)       = performance(cond);
            fitted_model2(beta).avsamples(cond)     = fitted_avsamples(cond);
            fitted_model2(beta).ll(cond)            = ll;

        end % end of conditions loop

    end % end of betas loop

    allsim_models2{1,cs} = simulated_model2;
    allfit_models2{1,cs} = fitted_model2;
end % end of Cs loop


%% COMPARE & SELECT MODELS %%



%% COMPUTE AQ DIFFERENCES %%

for sub = 1:nsubs

    % exctract this sub data
    sub_model1          = all_model1{1,sub};
    sub_model2          = all_model2{1,sub};

    % run the function
    AQdiffs     = computeAQdifference(sub_model2, condtrials, sub);

    all_AQdiffs{1,sub} = AQdiffs;


end % end of subjects loop


%%  PREP EEG DATA %%

% add eegpath
wd              = pwd;
eeg_path        = fullfile(wd, 'cropped'); addpath(genpath(eeg_path));


% run this first for ERPs and then for TFRs 
for sub = 1:nsubs 

    % load subject-specific cropped MEEG file
    sub_eeg         = load(fullfile(eeg_path, sprintf('erpcropped_data_sub_%02d.mat', sub)));
    sub_drawinfo    = allsub_drawinfo{1,sub};
    sub_AQdiffs     = all_AQdiffs{1,sub};
    sub_cond        = cond_data{sub,1};

end 

%% PLOT STUFF %%





if fparams == 2
    % extract free parameter
    param(1)    = R.initsample;
    param(2)    = R.initbeta;

else
    if R.modelnum == 1
        param       = R.initbeta;
    else
        param       = R.initsample;
    end
end




if R.freeparams == 2
    R.sample    = param(1);
    beta        = param(2); % softmax beta parameter is one of the free params
else
    if R.modelnum == 1 % if model is beta
        R.sample    = R.initsample;
        beta        = param; %
    else % if model is cost-to-sample
        R.sample    = param;
        beta        = R.initbeta; %
    end

end
