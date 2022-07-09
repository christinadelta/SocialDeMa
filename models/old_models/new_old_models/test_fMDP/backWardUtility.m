function Qsa = backWardUtility(sequence, draw, maxdraws, Cc, Cw, Cs, thisq)

utility = zeros(maxdraws, maxdraws+1);

ng = sum(sequence(1:draw));

for drawi = maxdraws : -1 : (draw + 1) %(backward - 10,9,8,7,6...
        
    [utility] = stateUtilityBeads(utility, drawi, draw, maxdraws, ng, Cc, Cw, Cs, thisq);
    
end
    
Qsa = actionValueBeads(utility, Cc, Cw, Cs, thisq, draw, ng, draw, maxdraws);
    
return