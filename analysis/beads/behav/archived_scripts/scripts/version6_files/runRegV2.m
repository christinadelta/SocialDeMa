function [fbetas, fcbetas, cpbetas, pbetas] = runRegV2(AQdiffs, allfrontal, allfc, allcp, allpar)

% RUN linear regressions using the fitlm function at the subject level for:
% - all frontal, left frontal, right frontal
% - all cf, left cf, right cf
% - all cp, left cp, right cp
% - all parietal, left parietal, right parietal

% OUTPUT: cells with subject reggression beta values to be used for
% one-sample t-tests

% this is only to test that both regression fucntions produce the same
% results. I mainly run runReg.m 

% -------------------
%%  extract arrays

% loop over laterality
for l = 1:3
    
    % run regressions for frontal
    fmodel      = fitlm(AQdiffs, allfrontal{l});
    fbetas{l}   = table2array(fmodel.Coefficients(2,1));
    
    % run regressions for frontocentral
    fcmodel     = fitlm(AQdiffs, allfc{l});
    fcbetas{l}  = table2array(fcmodel.Coefficients(2,1));
    
    % run regressions for centroparietal
    cpmodel     = fitlm(AQdiffs, allcp{l});
    cpbetas{l}  = table2array(cpmodel.Coefficients(2,1));
    
    % run regressions for frontal
    pmodel      = fitlm(AQdiffs, allpar{l});
    pbetas{l}   = table2array(pmodel.Coefficients(2,1));
    
end % end of laterality loop

end % end of function