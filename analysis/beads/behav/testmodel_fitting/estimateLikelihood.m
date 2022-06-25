function [ll] = estimateLikelihood(params, sequence, setData, fixedParams, findPick)

[ll, pickTrial, dQvec, ddec, aQvec] = estimateLikelihoodf(params, sequence, setData, fixedParams, findPick);

return