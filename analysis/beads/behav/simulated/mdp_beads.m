function [ll, Qsat, cprob] = mdp_beads(simR, drawSequence)

maxDraws            = size(drawSequence,2); 
k                   = simR.k;
beta                = simR.initbeta;
simR.sample         = simR.initsample;


% get number of trials/sequences 
ntrials         = max(size(drawSequence));
Qsat            = nan(ntrials,maxDraws,k);

% init log-likelihood
ll              = 0;

for trial = 1:ntrials

    Qsad            = zeros(maxDraws, k); % action values for this sequence will be stored here
    thisSequence    = drawSequence(trial,:);
    
    for draw = 1 : maxDraws
                   
        Qsad(draw, :) = backWardUtil(thisSequence, draw, maxDraws, simR);
        vVec          = Qsad(draw,:);

        % compute choices using softmax 
        cprobVals               = exp(beta*vVec)./sum(exp(beta*vVec)); % softmax to convert action values to probabiltiies of these actions
        cvec                    = cprobVals';
        cprob(trial,draw,:)     = cprobVals; 
        
    end

    % store this_trial Q values in a cell
    Qsat(trial,:,:) = Qsad;

    clear Qsad 
end % end of trials loop


return