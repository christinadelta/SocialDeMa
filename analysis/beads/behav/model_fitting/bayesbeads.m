function [mparams, lla, aQvec] = bayesbeads(sequence, choiceVec, trialinfo, alpha, Cw, cost_diff, Cs, sub)

% choice indecies (i don't think i need these)
b               = 1;
g               = 2;
s               = 3;

% extract trialinfo
q               = trialinfo.q;
urn             = trialinfo.urn;
cond            = trialinfo.cond;
accurate        = trialinfo.accurate;

% when running the func through prepro_beads.m comment this part 
% sequence        = this_sequence;
% choiceVec       = this_response;

% params          = Cw;
params          = cost_diff; 
fixedparams     = [alpha; q; Cs];
% findpick        = 0;
findpick        = 1;

% optimizing parameters for invidual subjects
options         = optimset('MaxFunEvals', 5000, 'TolFun', 0.001);
initialjitter   = [ 5; -5; 5; -5];

llamin          = Inf;

startparam      = params;

% run fminsearch
% lla = @(params) estimateLikelihood(params, sequence, choiceVec, fixedparams, findpick);
% mparams = fminsearch(lla, startparam, options);

[mparams, lla]  = fminsearch(@(params) estimateLikelihood(params, sequence, choiceVec, fixedparams, findpick),startparam, options);

if lla < llamin
    
    llamin      = lla;
    minparams   = mparams;
    
end

fprintf('ll %.3f\n', lla);

[ll, picktrial, dQvec, ddec, aQvec] = estimateLikelihoodf(minparams, sequence, choiceVec, fixedparams, findpick);



return