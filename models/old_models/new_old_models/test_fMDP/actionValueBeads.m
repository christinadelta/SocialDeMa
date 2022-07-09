function Qsa = actionValueBeads(utility, Cc, Cw, Cs, thisq, nd, ng, drawi, maxdraws)

pg = PG(thisq, nd, ng);

pb = 1 - pg;

QG = Cc*pg + Cw*pb;
QB = Cc*pb + Cw*pg;

if drawi < maxdraws

    QD = Cs + pb*((1-thisq)*utility(nd+1, ng+1+1) +   (thisq)*(utility(nd+1, ng+1))) + ...
                    pg*(  (thisq)*utility(nd+1, ng+1+1) + (1-thisq)*(utility(nd+1, ng+1)));
                
else
    
    QD = 0;
    
end

Qsa = [QG; QB; QD];

return