function prob_y = posteriorPredictive(x, postProb)

tvar = (1 + postProb.kappa)*postProb.sig/postProb.kappa;

sy = (x - postProb.mu)./sqrt(tvar);

prob_y = tpdf(sy, postProb.nu);

return