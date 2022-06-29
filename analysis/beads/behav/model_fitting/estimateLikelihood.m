function [ll] = estimateLikelihood(params, sequence, setData, fixedParams, findPick, urntype)

[ll, pickTrial, dQvec, ddec, aQvec, all_ll] = estimateLikelihoodf(params, sequence, setData, fixedParams, findPick, urntype);

return