function [mparams, lla, aQvec] = bayesbeads(sequence, choiceVec, info, alpha, Cw,Cc, cost_diff, Cs, cond, sub)

% extract trialinfo
thisq           = info.p;

% % % when running the func through prepro_beads.m comment this part 
% sequence        = thisub_seq;
% choiceVec       = thisub_choices;

% params          = Cw;
% params          = cost_diff; 
% fixedParams     = [alpha; thisq; cost_diff; cond];
fixedParams     = [alpha; thisq; Cw; Cc; cost_diff; cond];
params          = Cs;
findPick        = 1;

urntype         = info.urntypes;

options         = optimset('MaxFunEvals', 5000, 'TolFun', 0.001);
% initialJitter   = [ 5; -5; 5; -5];
% initialJitter   = [5   0.10;
%                  -5   0.10;
%                   5  -0.10;
%                  -5  -0.10];


llaMin          = Inf;

% for startValue = 1 : 4
    
%     initParam = params + initialJitter(1, :)'; 
%     startParam      = max(initParam);
    startParam = params;
    [mparams, lla] = fminsearch(@(params) estimateLikelihood(params, sequence, choiceVec, fixedParams, findPick, urntype),startParam, options);
    
    if lla < llaMin
        llaMin = lla;
        minParams = mparams;
    end

% end

fprintf('ll %.3f\n', lla);
fprintf('min params %.3f\n', minParams);

% ftxt = sprintf('subjectParams_%d_%d.mat', subject,types);
% save(ftxt, 'minParams', 'llaMin');

% [ll, pickTrial, dQvec, ddec, aQvec] = estimateLikelihoodf(minParams, sequence, choiceVec, fixedParams, findPick);
[ll, pickTrial, dQvec, ddec, aQvec] = estimateLikelihoodf(minParams, sequence, choiceVec, fixedParams, findPick, urntype);


return