function [v, d, Qvec] = Val(q, numDraws, numGreen, alpha, lseq, Cw, Cc, Cs)

% Computes action values Q for the 3 options (choose blue, choose
% green, draw again). 
% Maximum Q values determine the optimal action at each point in the
% sequence of draws. 

% compute probability that we  are drawing from green urn
pg          = PG(q, numDraws, numGreen); %  p(G|numDraws, numGreen)

% probability that we are drawing from predominantly blue urn
pb          = 1 - pg; % p(B|numDraws, numGreen) = 1 - p(G|numDraws, numGreen)

% determine the value of each available option. Here, only the cost of
% error is used. Cost of reward? Ok, Nick's 2011 paper explains why only Cw was
% used!!!!
QG          = Cw * pb + Cc * pg; % cost of choosing green 

QB          = Cw * pg + Cc * pb; % cost of choosing blue

% compute value for drawing again (QD)
if numDraws + 1 <= lseq
    
    % compute value of next state given that we draw a green bead
    val11   = vVal(q, numDraws+1, numGreen+1, alpha, lseq, Cw, Cc, Cs);
    
    % compute value of next state given that we draw a blue bead
    val10   = vVal(q, numDraws+1, numGreen, alpha, lseq, Cw, Cc, Cs);
    
    % action value is cost to sample plus expected value of future state
    QS      = Cs + pg * (val11 * q + val10 * (1 - q)) + pb * (val11 * (1 - q) + val10 * q);
    
    % arrange action values for the three options into a vector 
    Qvec    = [QG; QB; QS];
    
else % if it's the last sample then we only get action values for blue and green choices
    
    Qvec    = [QG; QB]; 
    
end

% softmax function to convert values to action probabilities -- mg(nd, ng)
d           = exp(alpha * Qvec)./sum(exp(alpha * Qvec));

% average value of this state if we take actions with probability d
v           = d'*Qvec;

end