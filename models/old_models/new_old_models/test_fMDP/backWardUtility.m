function Qsa = backWardUtility(drawSequence, draw, maxDraws, R)

utility = zeros(maxDraws, maxDraws+1);

ng = sum(drawSequence(1:draw));

for drawi = maxDraws : -1 : (draw + 1) %(backward - 10,9,8,7,6...
        
    [utility] = stateUtilityBeads(utility, drawi, draw, maxDraws, ng, R);
    
end
    
Qsa = actionValueBeads(utility, R, draw, ng, draw, maxDraws);
    
return