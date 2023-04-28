function [minParams, lla, Qsad, cprob, model_samples,model_urnchoice] = bayesbeads_b(thiscond_seqmat, cond_choices, R)

% RUN MODEL FITTING USING BRUNO'S VERSION

% this function calls fminsearch to fit free parameters:
%                       - model 1: free parameter 1 = Cost-sample
%                       - model 2: free parameter 1 = Cost-error
%                       - model 3: free parameter 1 = beta
%                       - model 4: free parameter 2 = Cost-sample, Cost-error
%                       - model 5: free parameter 2 = Cost-sample, beta

% POTENTIAL MODELS TO INCLUDE (will decide after making the above model fit and parameter recovery work):
%                       - model 6: free parameter 1 = discounting factor (gamma)
%                       - model 7: free parameter 3 = Cost-sample, beta, discounting factor (gamma)


% inputs: 1) m*k matrix (m= number of sequences, k=number of draws) --> i.e., 26x10
%         2) cell with choices (individual responses on every sequence)
%         every vector in the cell should be draws x options (or up to 10 x 3)
%         3) R structure with parameters and variables (e.g., cost-correct, cost-error,
%         cost-to-sample...)

% outputs: 1) mparams is the estimated optimal parameters (in our case
%          that is cost-to-sample and beta)
%          2) lla (minimum log likelihood)
%          3) Qsad cell with action values vectors.
%          4) cprob 3D matrix with choice probabilities for each action

% ------------------------------------------------------------------------------

% which parameter(s) to free?
model = R.model; % hich model?

if model == 1
    param       = R.initsample;
elseif model == 2
    % param     = R.initdiff;
    param(1)    = R.initerror;
    param(2)    = R.initreward;
elseif model == 3
    param(1)    = R.initsample;
    param(2)    = R.initerror;
    param(3)    = R.initreward;
elseif model == 4
    param       = R.initbeta;
elseif model == 5
    param(1)    = R.initsample;
    param(2)    = R.initbeta;
end

fitm            = 1; % when set to 1, fitmdp_beadsb runs to indicate that subjects' choice vector will be used 

% Define the optimization settings
options         = optimset('TolFun', 1e-6, 'MaxIter', 5000, 'Display', 'off');
% options         = optimset('Display','iter','MaxFunEvals', 5000, 'TolFun', 0.001,'PlotFcns', @optimplotfval);
% options         = optimset('MaxFunEvals', 5000, 'TolFun', 0.001); % the above is slow
    
llaMin          = Inf;
startParam      = param;

% Define a range of step sizes to try
% step_sizes      = [1, 0.1, 0.01, 0.001];

% for i = 1:length(step_sizes)

%     % Set the step size for this iteration
%     step_size = step_sizes(i);

    % Call the fminsearch function with the current step size
[mparams, lla]  = fminsearch(@(param) fitmdp_beadsb(param, R, thiscond_seqmat, cond_choices, fitm),startParam, options);

%     % Update the initial value for the next iteration
%     beta_init(i) = mparams

if lla < llaMin
    llaMin      = lla;
    minParams   = mparams;
end

%     % Adjust the step size for the next iteration
%     step_size = step_size * 0.1;
% end

% display fitted values on the command window
% if fparams == 2
%     fprintf('min ll: %.3f\n', lla);
%     fprintf('optimal cost-sample: %.3f\n ', minParams(1));
%     fprintf('optimal beta: %.3f\n ', minParams(2));
% else
%     fprintf('min ll: %.3f\n', lla);
%     % fprintf('optimal sample: %.3f\n ', minParams);
%     fprintf('optimal beta: %.3f\n ', minParams);
% end

fitm                                    = 0; % switch to zero to indicate that this will not use subject's choice vector.

[~, Qsad, cprob]                        = fitmdp_beads(minParams, R, thiscond_seqmat, cond_choices, fitm);

% compute model sampling rate based on choice probbilities via softmax 
N                                       = 1000; % number of iterations for computing model samples
[model_samples,model_urnchoice]         = computeModelSamples(cprob,N);


return