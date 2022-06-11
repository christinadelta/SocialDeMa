function [mparams, lla, aQvec] = bayesbeads(sequence, choiceVec, subject, types, p, cost_diff)

% choice indecies
b               = 1;
g               = 2;
s               = 3;

params          = [Cw ];
fixedparams     = [alpha; thisq; Cs];
findpick        = 0;

% optimizing parameters for invidual subjects
options         = optimset('MaxFunEvals', 5000, 'TolFun', 0.001);
initialjitter   = [ 5; -5; 5; -5];

llamin          = Inf;

startparam      = params;

[mparams, lla]  = fminsearch(@(params) estimateLikelihood(params, sequence, choiceVec, fixedparams, findpick),startparam, options);

if lla < llamin
    
    llamin      = lla;
    minparams   = mparams;
    
end

fprintf('ll %.3f\n', lla);

[ll, picktrial, dQvec, ddec, aQvec] = estimateLikelihoodf(minparams, sequence, choiceVec, fixedparams, findpick);


return