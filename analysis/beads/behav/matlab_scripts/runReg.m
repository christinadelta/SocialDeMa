function [fbetas, fcbetas, cpbetas, pbetas] = runReg(new_AQs, f, fc, cp, par)

% RUN linear regressions at the subject level for:
% - all frontal, left frontal, right frontal
% - all cf, left cf, right cf
% - all cp, left cp, right cp
% - all parietal, left parietal, right parietal

% OUTPUT: structs with subject reggression beta values to be used for
% one-sample t-tests

% -------------------
%%  extract arrays

% extract frontal arrays
frontal = f.f;
lf      = f.fl;
rf      = f.fr;

% extract fc arrays
fcent   = fc.fc;
lfc     = fc.lfc;
rfc     = fc.rfc;

% extract cp arrays
cpar    = cp.cp;
lcp     = cp.lcp;
rcp     = cp.rcp;

% extract parietal arrays
parietal    = par.p;
lp          = par.lp;
rp          = par.rp;

%%  run regressions

% run frontal regressions and store the beta coefficients 
% fmdl = fitlm(new_AQs, frontal)
fb = regress(frontal, new_AQs);
lfb = regress(lf, new_AQs);
rfb = regress(rf, new_AQs);

% run frontocentral regressions and store the beta coefficients
fcb = regress(fcent, new_AQs);
lfcb = regress(lfc, new_AQs);
rfcb = regress(rfc, new_AQs);

% run centroparietal regressions and store the beta coefficients
cpb = regress(cpar, new_AQs);
lcpb = regress(lcp, new_AQs);
rcpb = regress(rcp, new_AQs);

% run parietal regressions and store the beta coefficients
pb = regress(parietal, new_AQs);
lpb = regress(lp, new_AQs);
rpb = regress(rp, new_AQs);

% store all beta values in their structures

fbetas.f = fb;
fbetas.lf = lfb;
fbetas.rf = rfb;

fcbetas.fc = fcb;
fcbetas.lfc = lfcb;
fcbetas.rfc = rfcb;

cpbetas.cp = cpb;
cpbetas.lcp = lcpb;
cpbetas.rcp = rcpb;

pbetas.p = pb;
pbetas.lp = lpb;
pbetas.rp = rpb;

end % end of function