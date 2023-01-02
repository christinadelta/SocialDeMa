function AQdiffs = computeAQdiff(cond_modelQs, condtrials)

% this function is part of the BEADS task analyses 
% it computes difference in AQ values for each subject and for each condition
% OUTPUT: AQ difference values

% this function is part of the BEADS task analyses 

% trialQ_diff     = nan(1,1);
% counter         = 1;

% compute difference in AQ values for each subject and for each condition
for i = 1:condtrials
    
    % extract this_trial Q values
    trialQs                     = cond_modelQs{1,i};
    counter                     = 1;
    
    for j = 1:size(trialQs,1)
        
        thisQ                   = trialQs(j,:);                 % extract current draw's values 
        thisQ_diff              = thisQ(3) - max(thisQ(1:2));   % compute difference
        trialQ_diff(counter,1)    = thisQ_diff; 
        if thisQ_diff < 0
            break
        end
        counter                 = counter + 1;
    end
    
    % if trialqs and trialQ_diff not same length, add zeros to trialQ_diff
    if size(trialQs,1) ~= size(trialQ_diff,1)
        
        thisd = size(trialQs,1) - size(trialQ_diff,1)
        trialQ_diff(end+1:end+thisd) = nan
    end
    
    % store trialQ_diff
    AQdiffs{1,i} = trialQ_diff;
     
    clear trialQs counter trialQ_diff thisQ thisQ_diff

end

end