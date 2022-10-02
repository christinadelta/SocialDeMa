function [mparams, lla, Qsat] = fitbeads(thiscond_seqmat, thiscond_choiceVes, R)

% what is the free parameter? 
params = R.sample;

fixedParams = [R.correct; R.error; R.thisq; R.alpha];



options         = optimset('MaxFunEvals', 5000, 'TolFun', 0.001);

llaMin          = Inf;


return