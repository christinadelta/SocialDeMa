%% PARAMETER RECOVERY FOR MODEL FOR MODEL 1 (BETA) %%

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
%simR.initsample         = -0.25; % cost-to-sample is fixed in beta model
simR.freeparams         = 1; 
simR.initbeta           = 3;

Cs_bounds               = [-2 0]; % maximum and minimum value of beta?
nbins                   = 8;
allCs                   = linspace(Cs_bounds(1), Cs_bounds(2), nbins+1);

nReps                   = length(allCs);

% RUN MODEL 1 -- BETA
for rep = 1:nReps

    simR.initsample       = allCs(1,rep); % this rep beta value?
    
    % simulate dataset of 52 sequences/trials (26 easy and 26 difficult ones)
    simoutput           = simBeadsData(simvars,simR);

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
        simulated_model1(rep).beta(cond)        = simR.initbeta;
        simulated_model1(rep).cs(cond)          = simR.initsample;
        simulated_model1(rep).correct(cond)     = simoutput.performance(cond);

        % store fitted parameters and samples
        fitted_model1(rep).samples(:,cond)     = fitted_samples;
        fitted_model1(rep).Cs(cond)            = minParams(1);
        fitted_model1(rep).correct(cond)       = performance(cond);
        fitted_model1(rep).avsamples(cond)     = fitted_avsamples(cond);
        fitted_model1(rep).ll(cond)            = ll;
        
    end 

end % end of model 1 repetitions loop

clear simR fitted_samples fitted_urnchoice fitted_choice