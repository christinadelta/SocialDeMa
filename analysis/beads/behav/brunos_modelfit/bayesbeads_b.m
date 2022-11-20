function [minParams, ll, Qsad, cprob] = bayesbeads_b(thiscond_seqmat, thiscond_choiceVec, R)

% RUN MODEL FITTING USING BRUNO'S VERSION

% this function calls fminsearch to fit the free parameter cost-to-sample or Cs)
% inputs: 1) m*k matrix (m= number of sequences, k=number of draws) --> i.e., 26x10
%         2) cell with choices (individual responses on every sequence)
%         every vector in the cell should be draws x options (or up to 10 x 3)
%         3) R structure with parameters (e.g., cost-correct, cost-error,
%         cost-to-sample...)

% outputs: 1) mparams is the estimated fitted (free) parameters (in our case
%          that is cost-to-sample and beta)
%          2) lla (minimum log likelihood)
%          3) cell with action values vectors.

% ------------------------------------------------------------------------------

% extract free parameter
param(1) = R.initsample;
param(2) = R.initbeta;

% options         = optimset('Display','iter','MaxFunEvals', 5000, 'TolFun', 0.001,'PlotFcns', @optimplotfval);
options         = optimset('MaxFunEvals', 5000, 'TolFun', 0.001); % the above is slow
    
llaMin          = Inf;
startParam      = param;

[mparams, lla] = fminsearch(@(param) fitmdp_beadsb(param, R, thiscond_seqmat, thiscond_choiceVec),startParam, options);

if lla < llaMin
    llaMin = lla;
    minParams = mparams;
end

fprintf('min ll: %.3f\n', lla);
fprintf('optimal cost-sample: %.3f\n ', minParams(1));
fprintf('optimal beta: %.3f\n ', minParams(2));


[ll, Qsad, cprob] = fitmdp_beads(minParams, R, thiscond_seqmat, thiscond_choiceVec);


return