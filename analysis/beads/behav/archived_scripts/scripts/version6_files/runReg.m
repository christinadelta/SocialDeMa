function [fbetas, fcbetas, cpbetas, pbetas] = runReg(AQdiffs, allfrontal, allfc, allcp, allpar)

% RUN linear regressions using the regress function at the subject level for:
% - all frontal, left frontal, right frontal
% - all cf, left cf, right cf
% - all cp, left cp, right cp
% - all parietal, left parietal, right parietal

% OUTPUT: cells with subject reggression beta values to be used for
% one-sample t-tests

% -------------------
%%  extract arrays

% first compute an array of ones to used as the constant term (intercept)
% what is the length of the values?
AQs_len         = length(AQdiffs);
X               = ones(AQs_len,1); % the constant term
X(:,2)          = AQdiffs;

% loop over laterality
for l = 1:3
    
    % run regressions for frontal
    fmodel      = regress(allfrontal{l}, X);
    fbetas{l}   = fmodel(2);
    
    % run regressions for frontocentral
    fcmodel     = regress(allfc{l}, X);
    fcbetas{l}  = fcmodel(2);
    
    % run regressions for centroparietal
    cpmodel     = regress(allcp{l}, X);
    cpbetas{l}  = cpmodel(2); 
    
    % run regressions for frontal
    pmodel      = regress(allpar{l}, X);
    pbetas{l}   = pmodel(2);
    
end % end of laterality loop

end % end of function