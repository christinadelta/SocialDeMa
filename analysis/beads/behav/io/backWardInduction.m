function [reward, Qsat] = backWardInduction(thiscond_seqmat, R)

% how many trials per conditon?
Ntrials         = size(thiscond_seqmat,1);
maxDraws        = 10;

K               = 3; % choice options

reward          = zeros(Ntrials, 1);

Qsat            = zeros(Ntrials, maxDraws, K);

parfor trial = 1 : Ntrials
    
    Qsad            = zeros(maxDraws, K); % action values for this sequence will be stored here
    drawSequence    = thiscond_seqmat(trial,:);
    
%     % convert twos in the sequence to zeros (otherwise calculation of majority bead colour will be wrong - see backwardUtility.m)
%     seq_ones = find(drawSequence == 1);
%     seq_twos = find(drawSequence == 2);
%     drawSequence(seq_twos) = 0; 
%     
    for draw = 1 : maxDraws
                   
        Qsad(draw, :) = backWardUtility(drawSequence, draw, maxDraws, R);
        
    end
    
    % store this_trial Q values in a cell
    Qsat(trial,:,:) = Qsad;
    
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