function [NLLmat, ffX, ssX, ffSamples, ssSamples]  = reArrangeParams(mNLL,mSimX, mFitX,mFitSampels,mSimSamples,conditions,m)

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


end % end of subjects loop


end 