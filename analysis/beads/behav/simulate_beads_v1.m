% simulate data for the beads task and run ideal observer 

% craeted March 2024 
% to run as part of the beads analysis...

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

%% PERFORM SEARCH GRID 

% WHY PERFORM GRID SEARCH WITH SIMULATED DATA?
% To explore the parameter space and understand how changes in 
% parameters affect model predictions. This can help identify plausible 
% ranges for parameters and reveal how sensitive the model's outcomes are 
% to different parameter values.
% This preliminary exploration can inform you about the behavior of your 
% model under controlled conditions, where you know the "true" parameters 
% used to generate the data. It's a way to familiarize yourself with the 
% model's dynamics and to ensure that the grid search process and your 
% model implementation are working as expected.

% initial variables and parameters
init_sample     = -0.25;
init_beta       = 1;
bounds          = [-3 0;   % cost_sample range
                    0 5]; % beta range
bins            = [20 25];
nparam          = size(bins,2);

simalpha_bin    = init_sample/bounds(1,2)*bins(1);
simbeta_bin     = init_beta/bounds(2,2)*bins(2);
best_negLL      = inf(1,2); % initialize with a large number

% create the parameter space for the grid search
for iParam = 1:nparam
    range       = linspace(bounds(iParam,1),bounds(iParam,2),bins(iParam)+1);
    p{iParam}   = range(2:end); % stay just off the zero bounds
end

% estimate nll using different combinations of parameter values 
for cond = 1:conditions

    % Initialize NLL storage for the current condition
    nll(:,:,cond)       = inf(bins(1), bins(2));
    simvars.cond        = cond;
    mout                = runIdealObserver(simvars);
    %mout                = output{cond};
    sequence_matrix     = mout.simsequences;
    condition_choices   = mout.simchoicevec;
    simvars.urntype     = mout.simurns;
    simvars.thisq       = simvars.qvals(cond);

    for t = 1:bins(1) % cost-sample loop
        for tt = 1:bins(2) % beta loop

            param           = [p{1}(t), p{2}(tt)];
            negLL           = beads_lik_v1(param, simvars, sequence_matrix, condition_choices);
            
            % wtore NLL
            nll(t,tt,cond)  = negLL;
            
            % check for best NLL and update parameters
            if negLL < best_negLL(1,cond)
                best_negLL(1,cond)      = negLL;
                best_params(cond, :)    = param;
            end
        end % end of beta loop
    end % end of cost-sample function
end % end of conditions loop

%% plot NLL 

figure; 

% Loop over all conditions
for cond = 1:conditions
    
    % Select subplot
    subplot(1,2,cond);
    
    % Extract the NLL matrix for the current condition
    current_nll = nll(:,:,cond);
    
    % Plotting
    imagesc(p{1}, p{2}, current_nll');
    colorbar;
    xlabel('Cost-Sample');
    ylabel('Beta');
    fontsize(gcf,16,"points")
    title(sprintf('NEGLL: Probability Condition: %.1f', simvars.qvals(cond)));
    
    % Highlight the best parameters
    hold on;
    best_sample = best_params(cond, 1);
    best_beta   = best_params(cond, 2);
    plot(best_sample, best_beta, 'rp', 'MarkerSize', 16, 'MarkerFaceColor', 'r');
    hold off;
       
end

% Adjust subplot spacing if needed
sgtitle('NEGLL for All Conditions'); % Super title for the whole figure
fontsize(gcf,16,"points")

%% estimate marginal likelihoods of parameter values

totalPairs = conditions;

% bins for beta are 25, initialize the matrices to store marginal NLLs
aggregated_marginal_nll_sample  = zeros(20, totalPairs); 
aggregated_marginal_nll_beta    = zeros(25, totalPairs);

% loop over conditions to calculate marginal likelihoods
for cond = 1:conditions 

    % convert NLL to likelihoods for marginalization
    likelihoods = exp(-nll(:,:,cond));
    
    % marginalize likelihoods over parameters
    marginal_likelihood_sample  = sum(likelihoods, 2); % Sum over beta
    marginal_likelihood_beta    = sum(likelihoods, 1); % Sum over samples
    
    % convert marginal likelihoods to NLL
    marginal_nll_sample         = -log(marginal_likelihood_sample);
    marginal_nll_beta           = -log(marginal_likelihood_beta);
    
    % correctly aggregate the marginal NLLs
    aggregated_marginal_nll_sample(:, cond) = marginal_nll_sample;
    aggregated_marginal_nll_beta(:, cond)   = marginal_nll_beta(:); % Ensure column vector

end % end of conditions loop

% plot marginal cost-sample
% Plot Marginal NLL for Alpha across all conditions and volatilities
figure('Name', 'Marginal NLL for cost-sample across conditions', 'NumberTitle', 'off');

for cond = 1:conditions

    subplot(1,2,cond);
        
    % Plot for the current condition/volatility pair
    plot(p{1}, aggregated_marginal_nll_sample(:, cond), 'LineWidth', 2);
    xlabel('Cost-Sample');
    ylabel('Marginal Negative LL');
    title(sprintf('Probability Condition: %.1f', simvars.qvals(cond)));
    fontsize(gcf,14,"points")
    grid on;

end % end of codntiions loop

% plot marginal beta
% Plot Marginal NLL for beta across all conditions and volatilities
figure('Name', 'Marginal NLL for beta values across conditions', 'NumberTitle', 'off');

for cond = 1:conditions

    subplot(1,2,cond);
        
    % Plot for the current condition/volatility pair
    plot(p{2}, aggregated_marginal_nll_beta(:, cond), 'LineWidth', 2);
    xlabel('beta');
    ylabel('Marginal Negative LL');
    title(sprintf('Probability Condition: %.1f', simvars.qvals(cond)));
    fontsize(gcf,14,"points")
    grid on;

end % end of codntiions loop

%% estimate expected values of the parameter space 

% Initialize variables for expected values
expected_sample             = zeros(1, conditions);
expected_beta               = zeros(1, conditions);

% Loop over all condition/volatility pairs
for cond = 1:conditions
    
    % Convert aggregated marginal NLLs to probabilities for alpha
    prob_sample             = exp(-aggregated_marginal_nll_sample(:, cond) - min(aggregated_marginal_nll_sample(:, cond)));
    prob_sample             = prob_sample / sum(prob_sample);
    
    % Convert aggregated marginal NLLs to probabilities for beta
    prob_beta               = exp(-aggregated_marginal_nll_beta(:, cond) - min(aggregated_marginal_nll_beta(:, cond)));
    prob_beta               = prob_beta / sum(prob_beta);
    
    % Compute expected values for alpha and beta for the current pair
    expected_sample(cond)   = sum(p{1}(:) .* prob_sample(:));
    expected_beta(cond)     = sum(p{2}(:) .* prob_beta(:));
end

%% fit model with the best marginal and expected values from grid search

% Initialize storage for best marginal parameters
best_marginal_sample        = zeros(1, conditions);
best_marginal_beta          = zeros(1, conditions);

% create struct to store variables for fitting
R.ntrials         = totaltrials;
R.maxDraws        = 10;
R.qvals           = [0.8 0.6];
R.costsample      = -0.25;
R.correct         = 10;
R.error           = -10;
R.difference      = -20;
R.beta            = 1;
R.conditions      = conditions;
R.contrials       = totaltrials / conditions;

% Loop over all conditions to find best marginal parameters
for cond = 1:conditions

    % Find index of the best marginal alpha and beta (minimum NLL)
    [~, idx_sample]         = min(aggregated_marginal_nll_sample(:, cond));
    [~, idx_beta]           = min(aggregated_marginal_nll_beta(:, cond));
    
    % Store the best marginal alpha and beta values
    best_marginal_sample(cond)  = p{1}(idx_sample);
    best_marginal_beta(cond)    = p{2}(idx_beta);
    thiscond_model              = output{1,cond};

    % Fit the model with the best marginal parameter values for the current pair
    best_params                     = [best_marginal_sample(cond) best_marginal_beta(cond)];
    R.thisq                         = R.qvals(cond);
    modeloutput                     = fit_POMDP_Beads_v1(R, thiscond_model,best_params);
    marginal_model_choices(:,cond)  = modeloutput.samples;
    clear modeloutput

     % Fit the model with the expected parameter values for the current pair
    expected_params                 = [expected_sample(cond) expected_beta(cond)];
    modeloutput                     = fit_POMDP_Beads_v1(R, thiscond_model,expected_params);
    expected_model_choices(:,cond)  = modeloutput.samples;
    clear modeloutput
end % end of condition loop

%% %% plot marginal vs expected model choices (against true model choices)

% define dificulty level
probabilities = {'0.8', '0.6'};
figure;
for cond = 1:conditions

    subplot(1,2,cond);
    hold on;
    
    % Plot for marginal probabilities
    % Adjust indexing for Ps to match the 3x2 structure
    [~,~,h(1)] = myScatter(output{1,cond}.simdraws, marginal_model_choices(:,cond), false, [0 0 1], 'x');
    
    % Plot for expected probabilities
    [~,~,h(2)] = myScatter(output{1,cond}.simdraws, expected_model_choices(:,cond), false, [1 0 0], 'o');
    
    % Add unity line
    h(3) = plot([1 10], [1 10], 'k:', 'linewidth', 2);
    legend(h(1:2), {'Maximum Likelihood', 'Expected Value'}, 'location', 'best'); legend boxoff;
    xlim([1 10]); ylim([1 10]);
    xlabel('model choices - simulated');
    ylabel('model choices - estimated');

    % Title for each subplot with condition and volatility
    title(sprintf('%.1f Probability Condition',simvars.qvals(cond)));
    fontsize(gcf, 15, "points");
    
    % Display correlation coefficients for each subplot
    corrCoeffSimExpected = corrcoef(output{1,cond}.simdraws, expected_model_choices(:,cond));
    disp(['Condition ', num2str(cond), ' - Correlation (Simulated vs. Expected): ', num2str(corrCoeffSimExpected(1,2))]);
    
    corrCoeffMarginalExpected = corrcoef(marginal_model_choices(:,cond), expected_model_choices(:,cond));
    disp(['Condition ', num2str(cond), ' - Correlation (Marginal vs. Expected): ', num2str(corrCoeffMarginalExpected(1,2))]);
    
end  % end of condtions loop 

%% perform parameter recovery 

% In parameter recovery, this time I am randomly generating starting values

% how many repetitions?
repetitions         = 1000;

% define ranges for alpha and beta for true parameter values 
sample_range        = [-3 0]; % range for alpha
beta_range          = [0 5]; % range for beta

% create struct to store variables for fitting
R.ntrials           = totaltrials;
R.maxDraws          = 10;
R.qvals             = [0.8 0.6];
R.costsample        = -0.25;
R.correct           = 10;
R.error             = -10;
R.difference        = -20;
R.conditions        = conditions;
R.contrials         = totaltrials / conditions;

for rep = 1:repetitions

    for cond = 1:conditions

        % Generate random true values for alpha and beta within the specified ranges
        true_cs     = sample_range(1) + (sample_range(2) - sample_range(1)) * rand();
        true_beta   = beta_range(1) + (beta_range(2) - beta_range(1)) * rand();
        sim_params  = [true_cs true_beta];
        R.cond      = cond;
        R.thisq     = R.qvals(cond);

        % simulate data and model
        simoutput                   = sim_POMDP_Beads_v1(R,sim_params);

        % store data
        simX_sample(rep,cond)       = sim_params(1);
        simX_beta(rep,cond)         = sim_params(2);
        sim_samples(:,rep,cond)     = simoutput.simdraws;
        sim_sequences{rep,cond}     = simoutput.simsequences;
        sim_choices{rep,cond}       = simoutput.simchoicevec;
        R.urntype                   = simoutput.simurns;

        % start with model fitting
        fit_cs                      = sample_range(1) + (sample_range(2) - sample_range(1)) * rand();
        fit_beta                    = beta_range(1) + (beta_range(2) - beta_range(1)) * rand();

        % define the objective function
        obFunc                      = @(x) beads_lik_v1([x(1), x(2)], R, simoutput.simsequences, simoutput.simchoicevec);

        options = optimoptions('fmincon', 'Algorithm', 'interior-point', ...
                       'MaxIterations', 1000, ...
                       'OptimalityTolerance', 1e-6, ...
                       'StepTolerance', 1e-10, ...
                       'Display', 'iter'); % Shows iteration information

        X0                          = [fit_cs fit_beta];
        LB                          = [-3 0];
        UB                          = [0 5];
        [Xfit, NegLL]               = fmincon(obFunc, X0, [], [], [], [], LB, UB, [], options);

        % store fitted values
        fitX_sample(rep,cond)       = Xfit(1);
        fitX_beta(rep,cond)         = Xfit(2);
        fit_NLL(rep,cond)           = NegLL;

        % run model using the optimal values
        fit_output                  = fit_POMDP_Beads_v1(R, simoutput,Xfit);
        fit_samples(:,rep,cond)     = fit_output.samples;

 
    end % end of conditions loop
end % end of repetitions loop


%% plot cost-sample parameter values 

qvals = [0.8 0.6];

% Setting up the figure
figure;

% loop through each condition
for i = 1:conditions
    
    subplot(1, 2, i); % Create subplot for each condition
    scatter(simX_sample(:,i), fitX_sample(:,i)); % Scatter plot of simulated vs. fitted
    xlabel('Simulated Cost-Sample Values');
    ylabel('Fitted Cost-Sample Values');
    title(sprintf('%.1f Probability Condition', qvals(i)));
    
    % Calculate and display Spearman correlation
    [rho, pval] = corr(simX_sample(:,i), fitX_sample(:,i), 'Type', 'Spearman');
    corrText    = sprintf('rho = %.2f, p = %.3f', rho, pval); % Spearman correlation coefficient and p-value
    text(min(simX_sample(:,i)), max(fitX_sample(:,i)), corrText, 'VerticalAlignment', 'top'); % Display correlation

    % Calculate quantiles
    sim_quantiles = quantile(simX_sample(:,i), [0.25, 0.5, 0.75]);
    fit_quantiles = quantile(fitX_sample(:,i), [0.25, 0.5, 0.75]);
    quantileText = sprintf('Quantiles (Sim/Fit)\n25%%: %.2f/%.2f\n50%%: %.2f/%.2f\n75%%: %.2f/%.2f', ...
                            sim_quantiles(1), fit_quantiles(1), sim_quantiles(2), fit_quantiles(2), sim_quantiles(3), fit_quantiles(3));
    text(min(simX_sample(:,i)), min(fitX_sample(:,i)), quantileText, 'VerticalAlignment', 'bottom'); % Display quantiles

    % Calculate and plot line of best fit
    coeffs = polyfit(simX_sample(:,i), fitX_sample(:,i), 1); % First-degree polynomial coefficients for line of best fit
    xFit = linspace(min(simX_sample(:,i)), max(simX_sample(:,i)), 100); % Generate x values
    yFit = polyval(coeffs, xFit); % Calculate y values
    hold on; % Keep the scatter plot
    plot(xFit, yFit, 'r-', 'LineWidth', 2); % Plot line of best fit
    hold off;
end

sgtitle('Correlation and Quantiles for Simulated and Fitted Cost-Sample Parameter Values'); % Super title for the figure

%% plot beta parameter values 

qvals = [0.8 0.6];

% Setting up the figure
figure;

% loop through each condition
for i = 1:conditions
    
    subplot(1, 2, i); % Create subplot for each condition
    scatter(simX_beta(:,i), fitX_beta(:,i)); % Scatter plot of simulated vs. fitted
    xlabel('Simulated \beta Values');
    ylabel('Fitted \beta Values');
    title(sprintf('%.1f Probability Condition', qvals(i)));
    
    % Calculate and display Spearman correlation
    [rho, pval] = corr(simX_beta(:,i), fitX_beta(:,i), 'Type', 'Spearman');
    corrText    = sprintf('rho = %.2f, p = %.3f', rho, pval); % Spearman correlation coefficient and p-value
    text(min(simX_beta(:,i)), max(fitX_beta(:,i)), corrText, 'VerticalAlignment', 'top'); % Display correlation

    % Calculate quantiles
    sim_quantiles = quantile(simX_beta(:,i), [0.25, 0.5, 0.75]);
    fit_quantiles = quantile(fitX_beta(:,i), [0.25, 0.5, 0.75]);
    quantileText = sprintf('Quantiles (Sim/Fit)\n25%%: %.2f/%.2f\n50%%: %.2f/%.2f\n75%%: %.2f/%.2f', ...
                            sim_quantiles(1), fit_quantiles(1), sim_quantiles(2), fit_quantiles(2), sim_quantiles(3), fit_quantiles(3));
    text(min(simX_beta(:,i)), min(fitX_beta(:,i)), quantileText, 'VerticalAlignment', 'bottom'); % Display quantiles

    % Calculate and plot line of best fit
    coeffs = polyfit(simX_beta(:,i), fitX_beta(:,i), 1); % First-degree polynomial coefficients for line of best fit
    xFit = linspace(min(simX_beta(:,i)), max(simX_beta(:,i)), 100); % Generate x values
    yFit = polyval(coeffs, xFit); % Calculate y values
    hold on; % Keep the scatter plot
    plot(xFit, yFit, 'r-', 'LineWidth', 2); % Plot line of best fit
    hold off;
end

sgtitle('Correlation and Quantiles for Simulated and Fitted Cost-Sample Parameter Values'); % Super title for the figure

%% plot samples correlations

qvals = [0.8 0.6];

% Setting up the figure
figure;

for i = 1:conditions 

   % extract conditon parameter values
   cond_sim_samples = sim_samples(:,:,i); 
   cond_fit_samples = fit_samples(:,:,i); 
   av_sim_samples = mean(cond_sim_samples,2);
   av_fit_samples = mean(cond_fit_samples,2);

    subplot(1, 2, i); % Create subplot for each condition
    scatter(av_sim_samples, av_fit_samples); % Scatter plot of simulated vs. fitted
    xlabel('Simulated \beta Values');
    ylabel('Fitted \beta Values');
    title(sprintf('%.1f Probability Condition', qvals(i)));
    
    % Calculate and display Spearman correlation
    [rho, pval] = corr(av_sim_samples, av_fit_samples, 'Type', 'Spearman');
    corrText    = sprintf('rho = %.2f, p = %.3f', rho, pval); % Spearman correlation coefficient and p-value
    text(min(av_sim_samples), max(av_fit_samples), corrText, 'VerticalAlignment', 'top'); % Display correlation

    % Calculate quantiles
    sim_quantiles = quantile(av_sim_samples, [0.25, 0.5, 0.75]);
    fit_quantiles = quantile(av_fit_samples, [0.25, 0.5, 0.75]);
    quantileText = sprintf('Quantiles (Sim/Fit)\n25%%: %.2f/%.2f\n50%%: %.2f/%.2f\n75%%: %.2f/%.2f', ...
                            sim_quantiles(1), fit_quantiles(1), sim_quantiles(2), fit_quantiles(2), sim_quantiles(3), fit_quantiles(3));
    text(min(av_sim_samples), min(av_fit_samples), quantileText, 'VerticalAlignment', 'bottom'); % Display quantiles

    % Calculate and plot line of best fit
    coeffs = polyfit(av_sim_samples, av_fit_samples, 1); % First-degree polynomial coefficients for line of best fit
    xFit = linspace(min(av_sim_samples), max(av_sim_samples), 100); % Generate x
    yFit = polyval(coeffs, xFit); % Calculate y values
    hold on; % Keep the scatter plot
    plot(xFit, yFit, 'r-', 'LineWidth', 2); % Plot line of best fit
    hold off;

end % end of conditions loop 

sgtitle('Correlation between Simulated and Fitted Model Samples'); % Super title for the figure

%% perform parameter recovery for beta model

% In parameter recovery, this time I am randomly generating starting values

% how many repetitions?
repetitions         = 100;

% define ranges for alpha and beta for true parameter values 
% sample_range        = [-3 0]; % range for alpha
beta_range          = [0 5]; % range for beta

% create struct to store variables for fitting
R.ntrials           = totaltrials;
R.maxDraws          = 10;
R.qvals             = [0.8 0.6];
R.costsample        = -0.25;
R.correct           = 10;
R.error             = -10;
R.difference        = -20;
R.conditions        = conditions;
R.contrials         = totaltrials / conditions;

for rep = 1:repetitions

    for cond = 1:conditions

        % Generate random true values for alpha and beta within the specified ranges
        % true_cs     = sample_range(1) + (sample_range(2) - sample_range(1)) * rand();
        true_cs     = -0.25;
        true_beta   = beta_range(1) + (beta_range(2) - beta_range(1)) * rand();
        sim_params  = [true_cs true_beta];
        R.cond      = cond;
        R.urntype   = simoutput.simurns;
        R.thisq     = R.qvals(cond);

        % simulate data and model
        simoutput                   = sim_POMDP_Beads_v1(R,sim_params);

        % store data
        % simX_sample(rep,cond)       = sim_params(1);
        simX_beta_model(rep,cond)   = sim_params(2);
        sim_samples_beta(:,rep,cond)         = simoutput.simdraws;
        sim_sequences_beta{rep,cond}         = simoutput.simsequences;
        sim_choices_beta{rep,cond}           = simoutput.simchoicevec;

        % start with model fitting
        % fit_cs                      = sample_range(1) + (sample_range(2) - sample_range(1)) * rand();
        fit_beta                    = beta_range(1) + (beta_range(2) - beta_range(1)) * rand();

        % define the objective function
        obFunc                      = @(x) beads_lik_v1(x(1), R, simoutput.simsequences, simoutput.simchoicevec);
    
        X0                          = fit_beta;
        LB                          = 0;
        UB                          = 5;

        options                     = optimoptions('fmincon', 'Algorithm', 'interior-point', ...
                                       'MaxIterations', 1000, ...
                                       'OptimalityTolerance', 1e-6, ...
                                       'StepTolerance', 1e-10, ...
                                       'Display', 'iter'); % Shows iteration information

        [Xfit, NegLL]               = fmincon(obFunc, X0, [], [], [], [], LB, UB, [], options);

        % store fitted values
        % fitX_sample(rep,cond)       = Xfit(1);
        fitX_beta_model(rep,cond)         = Xfit(1);
        fit_NLL_beta(rep,cond)       = NegLL;

        % run model using the optimal values
        fit_output                  = fit_POMDP_Beads_v1(R, thiscond_model,Xfit);
        fit_samples_beta(:,rep,cond)     = fit_output.samples;

 
    end % end of conditions loop
end % end of repetitions loop


%% plot beta parameter values 

qvals = [0.8 0.6];

% Setting up the figure
figure;

% loop through each condition
for i = 1:conditions
    
    subplot(1, 2, i); % Create subplot for each condition
    scatter(simX_beta_model(:,i), fitX_beta_model(:,i)); % Scatter plot of simulated vs. fitted
    xlabel('Simulated \beta Values');
    ylabel('Fitted \beta Values');
    title(sprintf('%.1f Probability Condition', qvals(i)));
    
    % Calculate and display Spearman correlation
    [rho, pval] = corr(simX_beta_model(:,i), fitX_beta_model(:,i), 'Type', 'Spearman');
    corrText    = sprintf('rho = %.2f, p = %.3f', rho, pval); % Spearman correlation coefficient and p-value
    text(min(simX_beta_model(:,i)), max(fitX_beta_model(:,i)), corrText, 'VerticalAlignment', 'top'); % Display correlation

    % Calculate quantiles
    sim_quantiles = quantile(simX_beta_model(:,i), [0.25, 0.5, 0.75]);
    fit_quantiles = quantile(fitX_beta_model(:,i), [0.25, 0.5, 0.75]);
    quantileText = sprintf('Quantiles (Sim/Fit)\n25%%: %.2f/%.2f\n50%%: %.2f/%.2f\n75%%: %.2f/%.2f', ...
                            sim_quantiles(1), fit_quantiles(1), sim_quantiles(2), fit_quantiles(2), sim_quantiles(3), fit_quantiles(3));
    text(min(simX_beta_model(:,i)), min(fitX_beta_model(:,i)), quantileText, 'VerticalAlignment', 'bottom'); % Display quantiles

    % Calculate and plot line of best fit
    coeffs = polyfit(simX_beta_model(:,i), fitX_beta_model(:,i), 1); % First-degree polynomial coefficients for line of best fit
    xFit = linspace(min(simX_beta_model(:,i)), max(simX_beta_model(:,i)), 100); % Generate x values
    yFit = polyval(coeffs, xFit); % Calculate y values
    hold on; % Keep the scatter plot
    plot(xFit, yFit, 'r-', 'LineWidth', 2); % Plot line of best fit
    hold off;
end

sgtitle('Correlation and Quantiles for Simulated and Fitted Beta Parameter Values'); % Super title for the figure

%% plot samples correlations

qvals = [0.8 0.6];

% Setting up the figure
figure;

for i = 1:conditions 

   % extract conditon parameter values
   cond_sim_samples = sim_samples(:,:,i); 
   cond_fit_samples = fit_samples_beta(:,:,i); 
   av_sim_samples = mean(cond_sim_samples,2);
   av_fit_samples = mean(cond_fit_samples,2);

    subplot(1, 2, i); % Create subplot for each condition
    scatter(av_sim_samples, av_fit_samples); % Scatter plot of simulated vs. fitted
    xlabel('Simulated \beta Values');
    ylabel('Fitted \beta Values');
    title(sprintf('%.1f Probability Condition', qvals(i)));
    
    % Calculate and display Spearman correlation
    [rho, pval] = corr(av_sim_samples, av_fit_samples, 'Type', 'Spearman');
    corrText    = sprintf('rho = %.2f, p = %.3f', rho, pval); % Spearman correlation coefficient and p-value
    text(min(av_sim_samples), max(av_fit_samples), corrText, 'VerticalAlignment', 'top'); % Display correlation

    % Calculate quantiles
    sim_quantiles = quantile(av_sim_samples, [0.25, 0.5, 0.75]);
    fit_quantiles = quantile(av_fit_samples, [0.25, 0.5, 0.75]);
    quantileText = sprintf('Quantiles (Sim/Fit)\n25%%: %.2f/%.2f\n50%%: %.2f/%.2f\n75%%: %.2f/%.2f', ...
                            sim_quantiles(1), fit_quantiles(1), sim_quantiles(2), fit_quantiles(2), sim_quantiles(3), fit_quantiles(3));
    text(min(av_sim_samples), min(av_fit_samples), quantileText, 'VerticalAlignment', 'bottom'); % Display quantiles

    % Calculate and plot line of best fit
    coeffs = polyfit(av_sim_samples, av_fit_samples, 1); % First-degree polynomial coefficients for line of best fit
    xFit = linspace(min(av_sim_samples), max(av_sim_samples), 100); % Generate x
    yFit = polyval(coeffs, xFit); % Calculate y values
    hold on; % Keep the scatter plot
    plot(xFit, yFit, 'r-', 'LineWidth', 2); % Plot line of best fit
    hold off;

end % end of conditions loop 

sgtitle('Correlation between Simulated and Fitted Model Samples'); % Super title for the figure

%% run Parameter recovery for beta + Cs using fminsearch  

% In parameter recovery, this time I am randomly generating starting values

% how many repetitions?
repetitions         = 100;

% define ranges for alpha and beta for true parameter values 
sample_range        = [-3 0]; % range for alpha
beta_range          = [0 5]; % range for beta

% create struct to store variables for fitting
R.ntrials           = totaltrials;
R.maxDraws          = 10;
R.qvals             = [0.8 0.6];
R.costsample        = -0.25;
R.correct           = 10;
R.error             = -10;
R.difference        = -20;
R.conditions        = conditions;
R.contrials         = totaltrials / conditions;

for rep = 1:repetitions

    for cond = 1:conditions

        % Generate random true values for alpha and beta within the specified ranges
        true_cs     = sample_range(1) + (sample_range(2) - sample_range(1)) * rand();
        true_beta   = beta_range(1) + (beta_range(2) - beta_range(1)) * rand();
        sim_params  = [true_cs true_beta];
        R.cond      = cond;
        R.thisq     = R.qvals(cond);

        % simulate data and model
        simoutput                       = sim_POMDP_Beads_v1(R,sim_params);

        % store data
        simX_sample_v2(rep,cond)        = sim_params(1);
        simX_beta_v2(rep,cond)          = sim_params(2);
        sim_samples_v2(:,rep,cond)      = simoutput.simdraws;
        sim_sequences{rep,cond}         = simoutput.simsequences;
        sim_choices{rep,cond}           = simoutput.simchoicevec;
        R.urntype                       = simoutput.simurns;

        % start with model fitting
        fit_cs                          = sample_range(1) + (sample_range(2) - sample_range(1)) * rand();
        fit_beta                        = beta_range(1) + (beta_range(2) - beta_range(1)) * rand();
        starting_points                 = [fit_cs fit_beta];
        
        % define the objective function
        obFunc                          = @(x) beads_lik_v1([x(1), x(2)], R, simoutput.simsequences, simoutput.simchoicevec);

        options                         = optimset('TolFun', 1e-6, 'MaxIter', 5000, 'Display', 'off');

        
        [Xfit, NegLL]                   = fminsearch(obFunc,starting_points,options);

        % store fitted values
        fitX_sample_v2(rep,cond)        = Xfit(1);
        fitX_beta_v2(rep,cond)          = Xfit(2);
        fit_NLL_v2(rep,cond)            = NegLL;

        % run model using the optimal values
        fit_output                      = fit_POMDP_Beads_v1(R, simoutput,Xfit);
        fit_samples_v2(:,rep,cond)      = fit_output.samples;

 
    end % end of conditions loop
end % end of repetitions loop

%% plot cost-sample parameter values 

qvals = [0.8 0.6];

% Setting up the figure
figure;

% loop through each condition
for i = 1:conditions
    
    subplot(1, 2, i); % Create subplot for each condition
    scatter(simX_sample_v2(:,i), fitX_sample_v2(:,i)); % Scatter plot of simulated vs. fitted
    xlabel('Simulated Cost-Sample Values');
    ylabel('Fitted Cost-Sample Values');
    title(sprintf('%.1f Probability Condition', qvals(i)));
    
    % Calculate and display Spearman correlation
    [rho, pval] = corr(simX_sample_v2(:,i), fitX_sample_v2(:,i), 'Type', 'Spearman');
    corrText    = sprintf('rho = %.2f, p = %.3f', rho, pval); % Spearman correlation coefficient and p-value
    text(min(simX_sample_v2(:,i)), max(fitX_sample_v2(:,i)), corrText, 'VerticalAlignment', 'top'); % Display correlation

    % Calculate quantiles
    sim_quantiles = quantile(simX_sample_v2(:,i), [0.25, 0.5, 0.75]);
    fit_quantiles = quantile(fitX_sample_v2(:,i), [0.25, 0.5, 0.75]);
    quantileText = sprintf('Quantiles (Sim/Fit)\n25%%: %.2f/%.2f\n50%%: %.2f/%.2f\n75%%: %.2f/%.2f', ...
                            sim_quantiles(1), fit_quantiles(1), sim_quantiles(2), fit_quantiles(2), sim_quantiles(3), fit_quantiles(3));
    text(min(simX_sample_v2(:,i)), min(fitX_sample_v2(:,i)), quantileText, 'VerticalAlignment', 'bottom'); % Display quantiles

    % Calculate and plot line of best fit
    coeffs = polyfit(simX_sample_v2(:,i), fitX_sample_v2(:,i), 1); % First-degree polynomial coefficients for line of best fit
    xFit = linspace(min(simX_sample_v2(:,i)), max(simX_sample_v2(:,i)), 100); % Generate x values
    yFit = polyval(coeffs, xFit); % Calculate y values
    hold on; % Keep the scatter plot
    plot(xFit, yFit, 'r-', 'LineWidth', 2); % Plot line of best fit
    hold off;
end

sgtitle('Correlation and Quantiles for Simulated and Fitted Cost-Sample Parameter Values'); % Super title for the figure

%% plot beta parameter values 

qvals = [0.8 0.6];

% Setting up the figure
figure;

% loop through each condition
for i = 1:conditions
    
    subplot(1, 2, i); % Create subplot for each condition
    scatter(simX_beta_v2(:,i), fitX_beta_v2(:,i)); % Scatter plot of simulated vs. fitted
    xlabel('Simulated \beta Values');
    ylabel('Fitted \beta Values');
    title(sprintf('%.1f Probability Condition', qvals(i)));
    
    % Calculate and display Spearman correlation
    [rho, pval] = corr(simX_beta_v2(:,i), fitX_beta_v2(:,i), 'Type', 'Spearman');
    corrText    = sprintf('rho = %.2f, p = %.3f', rho, pval); % Spearman correlation coefficient and p-value
    text(min(simX_beta_v2(:,i)), max(fitX_beta_v2(:,i)), corrText, 'VerticalAlignment', 'top'); % Display correlation

    % Calculate quantiles
    sim_quantiles = quantile(simX_beta_v2(:,i), [0.25, 0.5, 0.75]);
    fit_quantiles = quantile(fitX_beta_v2(:,i), [0.25, 0.5, 0.75]);
    quantileText = sprintf('Quantiles (Sim/Fit)\n25%%: %.2f/%.2f\n50%%: %.2f/%.2f\n75%%: %.2f/%.2f', ...
                            sim_quantiles(1), fit_quantiles(1), sim_quantiles(2), fit_quantiles(2), sim_quantiles(3), fit_quantiles(3));
    text(min(simX_beta_v2(:,i)), min(fitX_beta_v2(:,i)), quantileText, 'VerticalAlignment', 'bottom'); % Display quantiles

    % Calculate and plot line of best fit
    coeffs = polyfit(simX_beta_v2(:,i), fitX_beta_v2(:,i), 1); % First-degree polynomial coefficients for line of best fit
    xFit = linspace(min(simX_beta_v2(:,i)), max(simX_beta_v2(:,i)), 100); % Generate x values
    yFit = polyval(coeffs, xFit); % Calculate y values
    hold on; % Keep the scatter plot
    plot(xFit, yFit, 'r-', 'LineWidth', 2); % Plot line of best fit
    hold off;
end

sgtitle('Correlation and Quantiles for Simulated and Fitted Cost-Sample Parameter Values'); % Super title for the figure

%% plot samples correlations

qvals = [0.8 0.6];

% Setting up the figure
figure;

for i = 1:conditions 

   % extract conditon parameter values
   cond_sim_samples = sim_samples_v2(:,:,i); 
   cond_fit_samples = fit_samples_v2(:,:,i); 
   av_sim_samples = mean(cond_sim_samples,2);
   av_fit_samples = mean(cond_fit_samples,2);

    subplot(1, 2, i); % Create subplot for each condition
    scatter(av_sim_samples, av_fit_samples); % Scatter plot of simulated vs. fitted
    xlabel('Simulated \beta Values');
    ylabel('Fitted \beta Values');
    title(sprintf('%.1f Probability Condition', qvals(i)));
    
    % Calculate and display Spearman correlation
    [rho, pval] = corr(av_sim_samples, av_fit_samples, 'Type', 'Spearman');
    corrText    = sprintf('rho = %.2f, p = %.3f', rho, pval); % Spearman correlation coefficient and p-value
    text(min(av_sim_samples), max(av_fit_samples), corrText, 'VerticalAlignment', 'top'); % Display correlation

    % Calculate quantiles
    sim_quantiles = quantile(av_sim_samples, [0.25, 0.5, 0.75]);
    fit_quantiles = quantile(av_fit_samples, [0.25, 0.5, 0.75]);
    quantileText = sprintf('Quantiles (Sim/Fit)\n25%%: %.2f/%.2f\n50%%: %.2f/%.2f\n75%%: %.2f/%.2f', ...
                            sim_quantiles(1), fit_quantiles(1), sim_quantiles(2), fit_quantiles(2), sim_quantiles(3), fit_quantiles(3));
    text(min(av_sim_samples), min(av_fit_samples), quantileText, 'VerticalAlignment', 'bottom'); % Display quantiles

    % Calculate and plot line of best fit
    coeffs = polyfit(av_sim_samples, av_fit_samples, 1); % First-degree polynomial coefficients for line of best fit
    xFit = linspace(min(av_sim_samples), max(av_sim_samples), 100); % Generate x
    yFit = polyval(coeffs, xFit); % Calculate y values
    hold on; % Keep the scatter plot
    plot(xFit, yFit, 'r-', 'LineWidth', 2); % Plot line of best fit
    hold off;

end % end of conditions loop 

sgtitle('Correlation between Simulated and Fitted Model Samples'); % Super title for the figure

%% perform parameter recovery for beta model using fminsearch

% In parameter recovery, this time I am randomly generating starting values

% how many repetitions?
repetitions         = 1000;

% define ranges for alpha and beta for true parameter values 
% sample_range        = [-3 0]; % range for alpha
beta_range          = [0 5]; % range for beta

% create struct to store variables for fitting
R.ntrials           = totaltrials;
R.maxDraws          = 10;
R.qvals             = [0.8 0.6];
R.costsample        = -0.25;
R.correct           = 10;
R.error             = -10;
R.difference        = -20;
R.conditions        = conditions;
R.contrials         = totaltrials / conditions;

for rep = 1:repetitions

    for cond = 1:conditions

        % Generate random true values for alpha and beta within the specified ranges
        % true_cs     = sample_range(1) + (sample_range(2) - sample_range(1)) * rand();
        true_cs     = -0.25;
        true_beta   = beta_range(1) + (beta_range(2) - beta_range(1)) * rand();
        sim_params  = [true_cs true_beta];
        R.cond      = cond;
        R.urntype   = simoutput.simurns;
        R.thisq     = R.qvals(cond);


        % simulate data and model
        simoutput                   = sim_POMDP_Beads_v1(R,sim_params);

        % store data
        % simX_sample(rep,cond)       = sim_params(1);
        simX_beta_model_v2(rep,cond) = sim_params(2);
        sim_samples_beta_v2(:,rep,cond)         = simoutput.simdraws;
        sim_sequences{rep,cond}         = simoutput.simsequences;
        sim_choices{rep,cond}           = simoutput.simchoicevec;

        % start with model fitting
        % fit_cs                      = sample_range(1) + (sample_range(2) - sample_range(1)) * rand();
        fit_beta                    = beta_range(1) + (beta_range(2) - beta_range(1)) * rand();
        starting_points             = fit_beta;

        % define the objective function
        obFunc                      = @(x) beads_lik_v1(x(1), R, simoutput.simsequences, simoutput.simchoicevec);
    
        options                         = optimset('TolFun', 1e-6, 'MaxIter', 5000, 'Display', 'off');

        
        [Xfit, NegLL]                   = fminsearch(obFunc,starting_points,options);

        % store fitted values
        % fitX_sample(rep,cond)       = Xfit(1);
        fitX_beta_model_v2(rep,cond)    = Xfit(1);
        fit_NLL_beta_v2(rep,cond)       = NegLL;

        % run model using the optimal values
        fit_output                  = fit_POMDP_Beads_v1(R, thiscond_model,Xfit);
        fit_samples_beta_v2(:,rep,cond)     = fit_output.samples;

 
    end % end of conditions loop
end % end of repetitions loop

%% plot beta parameter values 

qvals = [0.8 0.6];

% Setting up the figure
figure;

% loop through each condition
for i = 1:conditions
    
    subplot(1, 2, i); % Create subplot for each condition
    scatter(simX_beta_model_v2(:,i), fitX_beta_model_v2(:,i)); % Scatter plot of simulated vs. fitted
    xlabel('Simulated \beta Values');
    ylabel('Fitted \beta Values');
    title(sprintf('%.1f Probability Condition', qvals(i)));
    
    % Calculate and display Spearman correlation
    [rho, pval] = corr(simX_beta_model_v2(:,i), fitX_beta_model_v2(:,i), 'Type', 'Spearman');
    corrText    = sprintf('rho = %.2f, p = %.3f', rho, pval); % Spearman correlation coefficient and p-value
    text(min(simX_beta_model_v2(:,i)), max(fitX_beta_model_v2(:,i)), corrText, 'VerticalAlignment', 'top'); % Display correlation

    % Calculate quantiles
    sim_quantiles = quantile(simX_beta_model_v2(:,i), [0.25, 0.5, 0.75]);
    fit_quantiles = quantile(fitX_beta_model_v2(:,i), [0.25, 0.5, 0.75]);
    quantileText = sprintf('Quantiles (Sim/Fit)\n25%%: %.2f/%.2f\n50%%: %.2f/%.2f\n75%%: %.2f/%.2f', ...
                            sim_quantiles(1), fit_quantiles(1), sim_quantiles(2), fit_quantiles(2), sim_quantiles(3), fit_quantiles(3));
    text(min(simX_beta_model_v2(:,i)), min(fitX_beta_model_v2(:,i)), quantileText, 'VerticalAlignment', 'bottom'); % Display quantiles

    % Calculate and plot line of best fit
    coeffs = polyfit(simX_beta_model_v2(:,i), fitX_beta_model_v2(:,i), 1); % First-degree polynomial coefficients for line of best fit
    xFit = linspace(min(simX_beta_model_v2(:,i)), max(simX_beta_model_v2(:,i)), 100); % Generate x values
    yFit = polyval(coeffs, xFit); % Calculate y values
    hold on; % Keep the scatter plot
    plot(xFit, yFit, 'r-', 'LineWidth', 2); % Plot line of best fit
    hold off;
end

sgtitle('Correlation and Quantiles for Simulated and Fitted Beta Parameter Values'); % Super title for the figure

%% plot samples correlations

qvals = [0.8 0.6];

% Setting up the figure
figure;

for i = 1:conditions 

   % extract conditon parameter values
   cond_sim_samples = sim_samples_beta_v2(:,:,i); 
   cond_fit_samples = fit_samples_beta_v2(:,:,i); 
   av_sim_samples = mean(cond_sim_samples,2);
   av_fit_samples = mean(cond_fit_samples,2);

    subplot(1, 2, i); % Create subplot for each condition
    scatter(av_sim_samples, av_fit_samples); % Scatter plot of simulated vs. fitted
    xlabel('Simulated \beta Values');
    ylabel('Fitted \beta Values');
    title(sprintf('%.1f Probability Condition', qvals(i)));
    
    % Calculate and display Spearman correlation
    [rho, pval] = corr(av_sim_samples, av_fit_samples, 'Type', 'Spearman');
    corrText    = sprintf('rho = %.2f, p = %.3f', rho, pval); % Spearman correlation coefficient and p-value
    text(min(av_sim_samples), max(av_fit_samples), corrText, 'VerticalAlignment', 'top'); % Display correlation

    % Calculate quantiles
    sim_quantiles = quantile(av_sim_samples, [0.25, 0.5, 0.75]);
    fit_quantiles = quantile(av_fit_samples, [0.25, 0.5, 0.75]);
    quantileText = sprintf('Quantiles (Sim/Fit)\n25%%: %.2f/%.2f\n50%%: %.2f/%.2f\n75%%: %.2f/%.2f', ...
                            sim_quantiles(1), fit_quantiles(1), sim_quantiles(2), fit_quantiles(2), sim_quantiles(3), fit_quantiles(3));
    text(min(av_sim_samples), min(av_fit_samples), quantileText, 'VerticalAlignment', 'bottom'); % Display quantiles

    % Calculate and plot line of best fit
    coeffs = polyfit(av_sim_samples, av_fit_samples, 1); % First-degree polynomial coefficients for line of best fit
    xFit = linspace(min(av_sim_samples), max(av_sim_samples), 100); % Generate x
    yFit = polyval(coeffs, xFit); % Calculate y values
    hold on; % Keep the scatter plot
    plot(xFit, yFit, 'r-', 'LineWidth', 2); % Plot line of best fit
    hold off;

end % end of conditions loop 

sgtitle('Correlation between Simulated and Fitted Model Samples'); % Super title for the figure

