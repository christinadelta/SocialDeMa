function Qsa = backWardUtility(drawSequence, draw, maxDraws, R)

utility = zeros(maxDraws, maxDraws+1); % 10x11 mat

% define number of green beads
numGreen = sum(drawSequence(1:draw));

for i = maxDraws : -1 : (draw + 1) % start backwards 10 to 1
    
    utility = stateUtilityBeads(utility, i, draw, maxDraws, numGreen, R);
    
end % end of for loop

Qsa = actionValueBeads(utility, R, draw, numGreen, draw, maxDraws);

end