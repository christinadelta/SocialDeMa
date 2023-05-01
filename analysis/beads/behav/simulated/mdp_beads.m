function [ll, Qsat, cprob] = mdp_beads(simR, drawSequence)

maxDraws            = size(drawSequence,2); 
k                   = simR.k;
beta                = simR.initbeta


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
        vVec  = Qsad(draw, :);

        % compute choices using softmax 
        cprob_blue              = exp(beta * vVec(1,1)) ./ sum(exp(beta * vVec(1,:)));
        cprob_green             = exp(beta * vVec(1,2)) ./ sum(exp(beta * vVec(1,:)));
        cprob_draw              = exp(beta * vVec(1,3)) ./ sum(exp(beta * vVec(1,:)));
        cVec                    = [cprob_blue cprob_green cprob_draw]';
        
        % store all choice probabilities 
        cprob(trial,draw,1)     = cprob_blue;
        cprob(trial,draw,2)     = cprob_green;
        cprob(trial,draw,3)     = cprob_draw;
        
    end

    % store this_trial Q values in a cell
    Qsat(trial,:,:) = Qsad;

end % end of trials loop


return