function [mparams, lla, aQvec] = bayesbeads(thiscond_seqmat, thiscond_choiceVec, R)

% RUN MODEL FITTING USING NICK'S VERSION

% this function calls fminsearch to fit the free parameter cost-to-sample or Cs)
% inputs: 1) m*k matrix (m= number of sequences, k=number of draws) --> i.e., 26x10
%         2) cell with choices (individual responses on every sequence)
%         every vector in the cell should be draws x options (or up to 10 x 3)
%         3) R structure with parameters (e.g., cost-correct, cost-error,
%         cost-to-sample...)

% outputs: 1) mparams is the estimated fitted (free) parameter (in our case
%          that is cost-to-sample)
%          2) lla (log likelihood)
%          3) cell with action values vectors.

% ------------------------------------------------------------------------------

% unpack free parameter
params          = R.sample;

% fixed parameters 
fixedParams     = [R.alpha; R.thisq; R.error; R.correct];
findPick        = 1;

options         = optimset('MaxFunEvals', 5000, 'TolFun', 0.001);  
llaMin          = Inf;

startParam      = params;

[mparams, lla]  = fminsearch(@(params) estimateLikelihood(params, thiscond_seqmat, thiscond_choiceVec, fixedParams, findPick),startParam, options);

if lla < llaMin
    llaMin      = lla;
    minParams   = mparams;
end

fprintf('ll %.3f\n', lla);
fprintf('min params %.3f\n', minParams);


[ll, pickTrial, dQvec, ddec, aQvec] = estimateLikelihoodf(minParams, thiscond_seqmat, thiscond_choiceVec, fixedParams, findPick);


return