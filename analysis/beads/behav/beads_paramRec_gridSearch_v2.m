% parameter recovery with Grid search 
% Beads Analysis using a POMDP

%% housekeeping commands

clear all
clc

%% set figure-docking as default 

set(0,'DefaultFigureWindowStyle','docked')

%% simulate some data and run ideal observer to get some responses

% define some initial variables for the simulations
totaltrials             = 52;
conditions              = 2;
simvars.ntrials         = totaltrials;
simvars.maxDraws        = 10;
simvars.qvals           = [0.8 0.6];
simvars.costsample      = -0.25;
simvars.correct         = 10;
simvars.error           = -10;
simvars.difference      = -20;
simvars.beta            = 1;
simvars.conditions      = conditions;
simvars.contrials       = totaltrials / conditions;

% loop over conditions and simulate data
for cond = 1:conditions 
    simvars.cond        = cond;
    output{cond}        = runIdealObserver(simvars);
end % end of conditions loop

%% PERFORM SEARCH GRID with fminsearch for beta+cs model

% initial variables and parameters
init_sample     = -0.25;
init_beta       = 1;
bounds          = [-3 0;  % cost_sample range
                    0 5]; % beta range
bins            = [6 6];
nparam          = size(bins,2);
reps            = 40; % 40 simulated subjects
allcs           = [-2 -1.5 -1 -0.5 -0.25 -0.05];
allbetas        = [0.05 0.5 1 2 3.5 5];

best_negLL      = inf(1,2); % initialize with a large number

% create the parameter space for the grid search
p{1} = allcs;
p{2} = allbetas;

for rep = 1:reps
    for cond = 1:conditions

        % init some condition-wise parameters
        simvars.cond        = cond;
        simvars.thisq       = simvars.qvals(cond);

        % loop over cs values 
        for t = 1:bins(1) % cost-sample loop
            for tt = 1:bins(2) % betas loop

                param       = [p{1}(t), p{2}(tt)];

                % simulate data and model
                simoutput   = sim_POMDP_Beads_v1(simvars,param);
                
                % store data
                simX_sample{cond}(t,tt,rep)     = param(1);
                simX_beta{cond}(t,tt,rep)       = param(2);
                sim_samples{cond,rep}(:,t,tt)   = simoutput.simdraws;
                sim_sequences{cond,rep}{t,tt}   = simoutput.simsequences;
                sim_choices{cond,rep}{t,tt}     = simoutput.simchoicevec;
                simvars.urntype                 = simoutput.simurns;

                % define the objective function
                obFunc                          = @(x) beads_lik_v1([x(1), x(2)], simvars, simoutput.simsequences, simoutput.simchoicevec);
                options                         = optimset('TolFun', 1e-6, 'MaxIter', 5000, 'Display', 'off');
                [Xfit, NegLL]                   = fminsearch(obFunc,param,options);

                % store fitted values 
                fitX_sample_v2{cond}(t,tt,rep)  = Xfit(1);
                fitX_beta_v2{cond}(t,tt,rep)    = Xfit(2);
                fit_NLL_v2{cond}(t,tt,rep)      = NegLL;
        
                % run model using the optimal values
                fit_output                      = fit_POMDP_Beads_v1(simvars, simoutput,Xfit);
                fit_samples_v2{cond,rep}(:,t,tt)= fit_output.samples;

            end % end of beta loop
        end % end of cs loop
    end % end of condtions loop
end % end of repetitions loop

%% plot sampling rates as bar plots

% bins, conditions, etc..
nBins_cs        = 6; % num of cost-sample bins
nBins_beta      = 6; % num of beta bins
conditions      = 2; % conditions
reps            = 40; % num of repetitions
qvals           = [0.8 0.6];

% define the color codes in a cell array
colours          = {'#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b'};

% loop through each condition
for cond = 1:conditions

    figure; % create a new figure for each condition
    sgtitle(sprintf('Sampling Rates with Error Bars for %.1f Probability Condition', qvals(cond))); % super title for the figure
    
    % loop through each beta value
    for betaIndex = 1:nBins_beta
        subplot(2, 3, betaIndex); % arrange subplots in 2 rows, 3 cols
        means   = zeros(1, nBins_cs); % to store mean sampling rates for each cost-sample
        sems    = zeros(1, nBins_cs); % to store SEM for each cost-sample

        % collect data for each cost-sample parameter at this beta
        for csIndex = 1:nBins_cs
            allSamples = []; % initialize to collect all samples for this cs-beta combo

            % collect samples from all repetitions
            for rep = 1:reps
                data        = fit_samples_v2{cond, rep}(:, csIndex, betaIndex); % Extract data
                allSamples  = [allSamples, data]; % concatenate horizontally
            end

            % calculate the mean sampling rate and SEM for this cs-beta combo
            means(csIndex)  = mean(allSamples, 'all');
            sems(csIndex)   = std(allSamples, 0, 'all') / sqrt(numel(allSamples));
        end

        % plot bars for mean sampling rates with the color for this beta value
        bar(means, 'FaceColor', colours{betaIndex});
        hold on;

        % add error bars for SEM
        errorbar(1:nBins_cs, means, sems, 'k', 'linestyle', 'none');
        
        % set plot properties
        title(sprintf('Beta = %.2f', p{2}(betaIndex)));
        xlabel('Cost-Sample Parameter');
        ylabel('Mean Sampling Rate');
        xticklabels(arrayfun(@(x) sprintf('%.2f', x), p{1}, 'UniformOutput', false));
        xlim([0 nBins_cs + 1]);
        ylim([0 10])
        box on;

        fontsize(gcf, 13, "points");
        
        hold off; % release plot for next iteration
    end
end

%% plot correlations of simulated vs estimated sampling rates 

% some params to consider
nBins_cs        = 6;
nBins_beta      = 6;
conditions      = 2; % Or however many conditions you have
reps            = 40; % Number of repetitions

% define the color codes in a cell array
colours          = {'#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b'};
qvals            = [0.8 0.6];

% Loop for each condition - assuming you might want to do this for each condition
for cond = 1:conditions
    figure;
    sgtitle(sprintf('Correlations for %.1f Probability Condition', qvals(cond)));

    % Create subplots for each combination of beta and cost-sample parameters
    for tt = 1:nBins_beta % beta
        for t = 1:nBins_cs % cost-sample

            % Linear index for subplot position
            subplotIdx = (tt - 1) * nBins_cs + t;
            subplot(nBins_beta, nBins_cs, subplotIdx);

            % Initialize arrays to collect all data points for current param combo
            allSimSamples = [];
            allFitSamples = [];

            % Collect data across repetitions
            for rep = 1:reps
                simData = sim_samples{cond, rep}(:, t, tt);
                fitData = fit_samples_v2{cond, rep}(:, t, tt);

                allSimSamples(:,rep) = simData;
                allFitSamples(:,rep) = fitData;
            end

            % get mean over trials for each repetition
            meanSimData = mean(allSimSamples,1)';
            meanFitData = mean(allFitSamples,1)';

            % plot scatter plot with averaged over trial samples (40
            % datapoints for simulated and 40 for fitted)
            scatter(meanSimData, meanFitData, 'filled', 'MarkerFaceColor', colours{tt});
            hold on;

            % calculate and plot Spearman correlation
            rho = corr(meanSimData, meanFitData, 'Type', 'Spearman'); % Use original, non-jittered data for correlation
            title(sprintf('Beta=%.2f, CS=%.2f\nSpearman ρ=%.2f', p{2}(tt), p{1}(t), rho));

            % add a line of best fit using original, non-jittered data
            coeffs = polyfit(meanSimData, meanFitData, 1); % Use original data for line of fit
            xFit = linspace(min(meanSimData), max(meanSimData), 100);
            yFit = polyval(coeffs, xFit);
            plot(xFit, yFit, 'Color', colours{tt}, 'LineWidth', 2);

            hold off;

            % Adjust axes, labels, etc.
            xlabel('Simulated');
            ylabel('Estimated');
            fontsize(gcf, 13, "points");
            axis tight;
        end
    end
end

%% plot MAE simulated and estimated parameter values 

% some params to consider
nBins_cs        = 6;
nBins_beta      = 6;
conditions      = 2; 
reps            = 40; % Number of repetitions

% define the color codes in a cell array
colours          = {'#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b'};
qvals            = [0.8 0.6];

% PLOT MAE FOR COST-SAMPLE
% loop for each condition
for cond = 1:conditions
    figure; % create a new figure for each condition
    sgtitle(sprintf('MAE of Fitted vs. Simulated Cost-Sample Parameters for %.1f Probability Condition', qvals(cond)));
    
    for tt = 1:nBins_beta % Loop over beta values
        for t = 1:nBins_cs % Loop over cost-sample values
            subplotIdx = (tt - 1) * nBins_cs + t;
            subplot(nBins_beta, nBins_cs, subplotIdx);
            
            % Extract the constant simulated value for the current t, tt combination
            simValue_sample     = simX_sample{cond}(t,tt,1); % Assuming the same across all repetitions
            
            % Extract fitted values for the current condition, t, and tt
            fittedValues_sample = squeeze(fitX_sample_v2{cond}(t,tt,:));
            
            % Calculate MAE for cost-sample and beta parameters
            mae_sample          = mean(abs(fittedValues_sample - simValue_sample));
            
            % Plotting the distribution of fitted values vs. the simulated constant value
            % For simplicity, this example plots only the cost-sample parameter
            histogram(fittedValues_sample, 'FaceColor', colours{tt});
            hold on;
            xline(simValue_sample, 'r-', 'LineWidth', 2); % Simulated value for cost-sample parameter
            hold off;
            
            title(sprintf('CS: %.2f, Beta: %.2f\nMAE Sample: %.2f', p{1}(t), p{2}(tt), mae_sample));
            xlabel('Fitted Values');
            ylabel('Frequency');

            fontsize(gcf, 13, "points");
            
            % Ensure the plots do not overlap
            axis tight;
        end
    end
end

% PLOT MAE FOR BETA
% loop for each condition
for cond = 1:conditions
    figure; % create a new figure for each condition
    sgtitle(sprintf('MAE of Fitted vs. Simulated Beta Parameters for %.1f Probability Condition', qvals(cond)));
    
    for tt = 1:nBins_beta % Loop over beta values
        for t = 1:nBins_cs % Loop over cost-sample values
            subplotIdx = (tt - 1) * nBins_cs + t;
            subplot(nBins_beta, nBins_cs, subplotIdx);
            
            % Extract the constant simulated value for the current t, tt combination         
            simValue_beta       = simX_beta{cond}(t,tt,1); % For beta parameter
            
            % Extract fitted values for the current condition, t, and tt
            fittedValues_beta   = squeeze(fitX_beta_v2{cond}(t,tt,:));
            
            % Calculate MAE for cost-sample and beta parameters
            mae_beta            = mean(abs(fittedValues_beta - simValue_beta));
            
            % Plotting the distribution of fitted values vs. the simulated constant value
            % For simplicity, this example plots only the cost-sample parameter
            histogram(fittedValues_beta, 'FaceColor', colours{tt});
            hold on;
            xline(simValue_beta, 'r-', 'LineWidth', 2); % Simulated value for cost-sample parameter
            hold off;
            
            title(sprintf('CS: %.2f, Beta: %.2f\nMAE Beta: %.2f', p{1}(t), p{2}(tt), mae_beta));
            xlabel('Fitted Values');
            ylabel('Frequency');
            fontsize(gcf, 13, "points");
            
            % Ensure the plots do not overlap
            axis tight;
        end
    end
end

%% PERFORM SEARCH GRID with fminsearch for beta model

% initial variables and parameters
init_sample         = -0.25;
init_beta           = 1;
bounds_beta         = [0 5]; % beta range
bins_beta           = 6;
nparam_beta         = size(bins,1);
reps                = 40; % 40 simulated subjects
best_negLL          = inf(1,1); % initialize with a large number

% create the parameter space for the grid search
p                   = linspace(bounds_beta(1,1),bounds_beta(1,2),bins_beta(1)+1);

for rep = 1:reps
    for cond = 1:conditions

        % init some condition-wise parameters
        simvars.cond        = cond;
        simvars.thisq       = simvars.qvals(cond);

        for tt = 1:bins_beta % betas loop

            param       = [init_sample, p(tt)];

            % simulate data and model
            simoutput   = sim_POMDP_Beads_v1(simvars,param);
            
            % store data
            simX_sample{cond}(tt,rep)     = param(1);
            simX_beta{cond}(tt,rep)       = param(2);
            sim_samples{cond,rep}(:,tt)   = simoutput.simdraws;
            sim_sequences{cond,rep}{tt}   = simoutput.simsequences;
            sim_choices{cond,rep}{tt}     = simoutput.simchoicevec;
            simvars.urntype               = simoutput.simurns;

            % define the objective function
            starting_points                     = p(tt);

            % define the objective function
            obFunc                              = @(x) beads_lik_v1(x(1), simvars, simoutput.simsequences, simoutput.simchoicevec);
            options                             = optimset('TolFun', 1e-6, 'MaxIter', 5000, 'Display', 'off');
            [Xfit, NegLL]                       = fminsearch(obFunc,starting_points,options);

            % store fitted values 
            fitX_beta_model_gs{cond}(tt,rep)    = Xfit(1);
            fit_NLL_beta_gs{cond}(tt,rep)       = NegLL;
    
            % run model using the optimal values
            fit_output                          = fit_POMDP_Beads_v1(simvars, simoutput,Xfit);
            fit_samples_beta_gs{cond,rep}(:,tt) = fit_output.samples;

        end % end of beta loop
    end % end of condtions loop
end % end of repetitions loop

