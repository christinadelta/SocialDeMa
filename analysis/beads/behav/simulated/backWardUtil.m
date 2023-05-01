function Qsa = backWardUtil(drawSequence, draw, maxDraws, R)

utility = zeros(maxDraws, maxDraws+1);

ng = sum(drawSequence(1:draw));

for drawi = maxDraws : -1 : (draw + 1)
        
    [utility] = stateUtilBeads(utility, drawi, draw, maxDraws, ng, R);
    
end
    
Qsa = compActionValues(utility, R, draw, ng, draw, maxDraws);
    
return