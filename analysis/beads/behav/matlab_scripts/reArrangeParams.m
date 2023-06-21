function [outmat, ffX, ssX] = reArrangeParams(allbetas, allCs, NLLcell, xfit, xsim, conditions)

% created in June 2023 
% The function takes as an input the cells from parameter recovery and
% re-arranges the data as a Cs-by-betas matrix (seperately for each
% condition

lenCs       = length(allCs);
lenBetas    = length(allbetas);
model       = 2; % re-arrangement is done for the CsBeta model

NLLdata     = NLLcell{1,model};
Xfitdata    = xfit{1, model};
Xsimdata    = xsim{1, model};

count       = 1;
%% Re-arrange NLLs

% loop over Cs 
for i = 1:lenCs

    d = NLLdata{1,i};
    
    % loop over conditions
    for j = 1:conditions 

        dc                  = d(:,j)';
        outmat{1,j}(i,:)    = dc;

    end % end of conditions loop
end % end of cost-sample 

%% Re-arrange params 

% loop over Cs 
for i = 1:lenCs

    sX = Xsimdata{1,i};
    fX = Xfitdata{1,i};

    % loop over betas
    for j = 1:lenBetas
        
        % loop over condition
        for k = 1:conditions

            ssX{1,k}(1,count) = sX.sample(j,k);
            ssX{1,k}(2,count) = sX.beta(j,k);

            ffX{1,k}(1,count) = fX.sample(j,k);
            ffX{1,k}(2,count) = fX.beta(j,k);

        end % end of condition loop

        count = count+1; % update counter

    end % end of betas loop

end % end of cs loop


end 