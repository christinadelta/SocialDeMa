function [reward, Qsat] = backWardInduction(Ntrials, maxDraws, allsequences, R)

K = 3;

reward = zeros(Ntrials, 1); % empty array 

Qsat = zeros(Ntrials, maxDraws, K); % Qsat 

parfor trial = 1 : Ntrials
    
    Qsad = zeros(maxDraws, 3); % this is like a choice vector?
    
    for draw = 1 : maxDraws
                   
        Qsad(draw, :) = backWardUtility(allsequences(trial, :), draw, maxDraws, R);
        
    end
    
    Qsat(trial, :, :) = Qsad;
    
    %%% randomize choice for symmetric values
    Qsac = Qsad + 0.000001*randn(maxDraws, K);
    
    Qsa1 = Qsac(:, 1) - Qsac(:, 3);
    Qsa2 = Qsac(:, 2) - Qsac(:, 3);
    
    choice1 = find(Qsa1 > 0);
    choice2 = find(Qsa2 > 0);
    
    if isempty(choice1)
        choice1(maxDraws+1) = 1;
    end
    
    if isempty(choice2)
        choice2(maxDraws+1) = 1;
    end
    
    if choice1(1) < choice2(1)
        reward(trial) = 1;
    else
        reward(trial) = 0;
    end
    
            
end
return