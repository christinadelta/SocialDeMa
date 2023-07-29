function refit_samples  = reArrangeParams(mdl_fitsamples)

% created in June 2023 
% The function takes as an input the cells from parameter recovery,
% and re-arranges the data for plotting

% output:
%           - fitsamples{1,cond}{1,1betas}(cs,reps)

% -------------------
% how many subjects?
betas = length(mdl_fitsamples);

% loop over betas 
for i = 1: betas

    thisb                           = mdl_fitsamples{1,i};

    for cond = 1:2

        tmp                         = thisb(:,:,cond);
        refit_samples{1,cond}{1,i}  = tmp;

    end % end of conditions loop

end % ned of betas loop


end % end of function