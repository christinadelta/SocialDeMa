function outmat = reArrangeParams(allbetas, allCs, incell, conditions)

% created in June 2023 
% The function takes as an input the cells from parameter recovery and
% re-arranges the data as a Cs-by-betas matrix (seperately for each
% condition

lenCs       = length(allCs);
lenBetas    = length(allbetas);
model       = 2; % re-arrangement is done for the CsBeta model

data        = incell{1,model};

% loop over Cs 
for i = 1:lenCs

    d = data{1,i};
    
    % loop over conditions
    for j = 1:conditions 

        dc                  = d(:,j)';
        outmat{1,j}(i,:)    = dc;

    end % end of conditions loop
end % end of cost-sample 

end 