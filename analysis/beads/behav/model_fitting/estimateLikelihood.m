function [ll] = estimateLikelihood(params, sequence, setdata, fixedparams, findpick)

% comment this out when running through prepro_beads.m 
setdata = choiceVec;

[ll, picktrial, dQvec, ddec, aQvec] = estimateLikelihoodf(params, sequence, setdata, fixedparams, findpick);


end