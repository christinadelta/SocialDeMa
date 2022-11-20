function AQdiffs = computeAQdiff(cond_modelQs, condtrials)

% this function is part of the BEADS task analyses 
% it computes difference in AQ values for each subject and for each condition
% OUTPUT: AQ difference values

% this function is part of the BEADS task analyses 

trialQ_diff     = nan(1,1);
counter         = 1;

% compute difference in AQ values for each subject and for each condition
for i = 1:condtrials
    
    % extract this_trial Q values
    trialQs                     = cond_modelQs{1,i};
    
    for j = 1:size(trialQs,1)
        
        thisQ                   = trialQs(j,:);                 % extract current draw's values 
        thisQ_diff              = thisQ(3) - max(thisQ(1:2));   % compute difference
        trialQ_diff(counter)    = thisQ_diff;                   
        counter                 = counter + 1;
    end
    
    clear trialQs

end

AQdiffs                         = trialQ_diff;
    

end