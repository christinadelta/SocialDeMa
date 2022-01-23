function [v, d, Qvec] = Val(q, numDraws, numGreen, alpha, sequenceL, Cw, Cs)

% compute probability that we  are drawing from green urn
pg          = PG(q, numDraws, numGreen);

% probability that we are drawing from predominantly blue urn
pb          = 1 - pg;

QG          = Cw * pb; % cost of choosing green 

QB          = Cw * pg; % cost of choosing blue

if numDraws + 1 <= sequenceL
    
    % compute value of next state given that we draw a green bead
    val11   = vVal(q, numDraws+1, numGreen+1, alpha, sequenceL, Cw, Cs);
    
    % compute value of next state given that we draw a blue bead
    val10   = vVal(q, numDraws+1, numGreen, alpha, sequenceL, Cw, Cs);
    
    % action value is cost to sample plus expected value of next state
    QS      = Cs + pg * (val11 * q + val10 * (1 - q)) + pb * (val11 * (1 - q) + val10 * q);
    
    % arrange action values into a vector 
    Qvec    = [QG; QB; QS];
    
else % if it's the last sample then we only get action values for blue and green choices
    
    Qvec    = [QG; QB]; 
    
end

% softmax function to convert values to action probabilities
d           = exp(alpha * Qvec)./sum(exp(alpha * Qvec));

% average value of this state if we take actions with probability d
v           = d'*Qvec;

end