function [ll] = estimateLikelihood_n(params, sequence, setData, fixedParams, findPick)

[ll, pickTrial, dQvec, ddec, aQvec] = estimateLikelihoodf_n(params, sequence, setData, fixedParams, findPick);

return