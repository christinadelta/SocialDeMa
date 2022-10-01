function [ll] = estimateLikelihood(params, sequence, setData, fixedParams, findPick, urntype)

[ll, pickTrial, dQvec, ddec, aQvec] = estimateLikelihoodf(params, sequence, setData, fixedParams, findPick, urntype);

return