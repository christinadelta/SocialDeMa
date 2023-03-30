function [minParams, lla, Qsad, cprob, model_samples,model_urnchoice] = bayesbeads_b(thiscond_seqmat, cond_choices, R)

% RUN MODEL FITTING USING BRUNO'S VERSION

% this function calls fminsearch to fit free parameters:
%                       - model 1: free parameter = beta
%                       - model 2: free parameter = beta, cost-to-sample

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

% how many free parameters?
fparams         = R.freeparams;

if fparams == 2
    % extract free parameter
    param(1)    = R.initsample;
    param(2)    = R.initbeta;

else
    param       = R.initbeta; 
end

fitm            = 1; % when set to 1, fitmdp_beadsb runs to indicate that subjects' choice vector will be used 

% options         = optimset('Display','iter','MaxFunEvals', 5000, 'TolFun', 0.001,'PlotFcns', @optimplotfval);
options         = optimset('MaxFunEvals', 5000, 'TolFun', 0.001); % the above is slow
    
llaMin          = Inf;
startParam      = param;

[mparams, lla]  = fminsearch(@(param) fitmdp_beadsb(param, R, thiscond_seqmat, cond_choices, fitm),startParam, options);

if lla < llaMin
    llaMin      = lla;
    minParams   = mparams;
end

if fparams == 2
    fprintf('min ll: %.3f\n', lla);
    fprintf('optimal cost-sample: %.3f\n ', minParams(1));
    fprintf('optimal beta: %.3f\n ', minParams(2));
else
    fprintf('min ll: %.3f\n', lla);
    fprintf('optimal beta: %.3f\n ', minParams);
end

fitm                                    = 0; % switch to zero to indicate that this will not use subject's choice vector.


[ll, Qsad, cprob]                       = fitmdp_beads(minParams, R, thiscond_seqmat, cond_choices, fitm);

% compute model sampling rate based on choice probbilities via softmax 
N                                       = 1000; % number of iterations for computing model samples
[model_samples,model_urnchoice]         = computeModelSamples(cprob,N);


return