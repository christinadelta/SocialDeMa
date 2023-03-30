function Qsa = backWardUtility_b(drawSequence, draw, maxDraws, R)

utility = zeros(maxDraws, maxDraws+1);

ng = sum(drawSequence(1:draw));

for drawi = maxDraws : -1 : (draw + 1)
        
    [utility] = stateUtilityBeads_b(utility, drawi, draw, maxDraws, ng, R);
    
end
    
Qsa = actionValueBeads_b(utility, R, draw, ng, draw, maxDraws);