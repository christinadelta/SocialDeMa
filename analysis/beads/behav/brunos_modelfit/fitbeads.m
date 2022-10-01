function [mparams, lla, Qsat] = fitbeads(thisub_choices, thisub_seq, info, R)

% what is the free parameter? 
params = R.sample;


options         = optimset('MaxFunEvals', 5000, 'TolFun', 0.001);

llaMin          = Inf;


return