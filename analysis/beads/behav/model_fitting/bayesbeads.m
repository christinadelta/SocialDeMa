function [mparams, lla, all_ll, aQvec] = bayesbeads(sequence, choiceVec, info, alpha, Cw, cost_diff, Cs, cond, sub)

% extract trialinfo
q               = info.p;

% % % when running the func through prepro_beads.m comment this part 
% sequence        = cond_sequence;
% choiceVec       = cond_choices;

% params          = Cw;
% params          = cost_diff; 
params          = Cs;
fixedParams     = [alpha; q; cost_diff; cond];
findPick        = 1;

urntype         = info.urntypes;

options         = optimset('MaxFunEvals', 5000, 'TolFun', 0.001);
initialJitter   = [ 5; -5; 5; -5];

llaMin          = Inf;
startParam      = params;

[mparams, lla] = fminsearch(@(params) estimateLikelihood(params, sequence, choiceVec, fixedParams, findPick, urntype),startParam, options);
    
if lla < llaMin
    llaMin = lla;
    minParams = mparams;
end

fprintf('ll %.3f\n', lla);

% ftxt = sprintf('subjectParams_%d_%d.mat', subject,types);
% save(ftxt, 'minParams', 'llaMin');

% [ll, pickTrial, dQvec, ddec, aQvec] = estimateLikelihoodf(minParams, sequence, choiceVec, fixedParams, findPick);
[ll, all_ll, pickTrial, dQvec, ddec, aQvec] = estimateLikelihoodf(minParams, sequence, choiceVec, fixedParams, findPick, urntype);


return