function [avNLL, avSimSample, avFitSample]  = reArrangeParams(mNLL,mSimX, mFitX,mFitSampels,mSimSamples,conditions,m)

% created in June 2023 
% The function takes as an input the cells from parameter recovery,
% avverages across iterations and re-arranges the data for visualisation

% how many subjects?
nsubs = size(mNLL,2);

% loop over participants 
for sub = 1:nsubs
    
    % extract sub parameters
    subNLL          = mNLL{1,sub};
    subSimX         = mSimX{1,sub};
    subFitX         = mFitX{1,sub};
    subFitSample    = mFitSampels{1,sub};
    subSimSample    = mSimSamples{1,sub};

    % loop over conditions 
    for cond = 1:conditions

        condNLL                 = subNLL(:,cond);
        avNLL(sub,cond)         = mean(condNLL);

        condSimSample           = subSimSample(:,cond);
        avSimSample(sub,cond)   = mean(condSimSample);

        condFitSample           = subFitSample(:,cond);
        avFitSample(sub,cond)   = mean(condFitSample);

        % now average the free params (TODO)



    end % end of conditions loop


end % end of subjects loop


end 