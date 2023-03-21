function AQdiffs = computeAQdiff_v2(sub_modelQs,condtrials,sub)

% this function is part of the BEADS task analyses 
% it computes difference in AQ values for each subject and for each condition
% and checks whether AQ values and participant draws are of the same length

% INPUT:     - moldel Q values
%            - number of condition trials
%            - subject number of draws

% OUTPUT: AQ difference values

% ----------------

trialQ_diff     = nan(1,1);
counter         = 1;
conds           = 2;

for cond = 1:conds
    
    % extrcat this_sub model Q values
    cond_Qs     = sub_modelQs{1,cond};
    
    % loop over trials now to compute differences in AQs
    for trl = 1:condtrials
        
        % extract this_trial Q values
        trialQs         = cond_Qs{1,trl};
        thisQ_len       = size(trialQs,1);
        
        % if this is sub 36 remove trial 14 of condition 2 because some of
        % the draws were not recorded by the biosemi device (for this trial
        % there are more draws than epochs)
        if sub == 36 
            if cond == 2 && trl == 14
            
                trialQs(1:5,:) = [];
                thisQ_len       = size(trialQs,1);
            end
        end
 
        % compute AQ difference for each draw
        for i = 1:thisQ_len
            thisQ                       = trialQs(i,:);                 % extract current draw's values 
            thisQ_diff                  = thisQ(3) - max(thisQ(1:2));   % compute difference
            trialQ_diff(counter,1)      = thisQ_diff; 
            counter                     = counter + 1;                  % update counter

        end % end of trial Qs loop

    end % end of condition trials
 
end % end of conditions loop

AQdiffs                                     = trialQ_diff; 

end