function [reward, Qsat] = backWardInduction(ntrials, maxDraws, drawSequence, R)

k       = 3;
reward  = zeros(ntrials,1); % create rewards array
Qsat    = zeros(ntrials, maxDraws, k); % why 3-dimensional?
 
for trial = 1:ntrials
    
    thisSequence    = drawSequence(trial,:);
    
    Qsad            = zeros(maxDraws,3); % is that related to behavioural choices? (choose green, choose blue, sample again)
    
    for draw = 1:maxDraws
        
        % run backwardutility function
        Qsad(draw, :) = backWardUtility(thisSequence, draw, maxDraws, R);
        
        
    end % end of draws loop 
end % end of trials loop

end