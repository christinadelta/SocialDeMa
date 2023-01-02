for i = 1:condtrials
    
    % extract this_trial Q values
    trialQs             = cond_modelQs{1,i};
    
    for j = 1:size(trialQs,1)
        
        thisQ           = trialQs(j,:); %
        thisQ_diff      = thisQ(3) - max(thisQ(1:2));
        trialQ_diff(j)  = thisQ_diff;
        
    end
    
    % store Q differences 
    AQdiffs{1,i}        = trialQ_diff;
    
    clear trialQ_diff
    
end