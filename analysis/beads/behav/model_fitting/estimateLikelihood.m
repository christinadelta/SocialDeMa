function [ll] = estimateLikelihood(params, sequence, setData, fixedparams, findpick)

[ll, picktrial, dQvec, ddec, aQvec] = estimateLikelihoodf(params, sequence, setdata, fixedparams, findpick);


end