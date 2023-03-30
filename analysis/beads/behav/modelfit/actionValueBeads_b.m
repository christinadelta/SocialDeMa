function Qsa = actionValueBeads_b(utility, R, nd, ng, drawi, maxDraws)

pg = PG_b(R.thisq, nd, ng);

pb = 1 - pg;

QG = R.correct*pg + R.error*pb;
QB = R.correct*pb + R.error*pg;

if drawi < maxDraws

    QD = R.sample + pb*((1-R.thisq)*utility(nd+1, ng+1+1) +   (R.thisq)*(utility(nd+1, ng+1))) + ...
                    pg*(  (R.thisq)*utility(nd+1, ng+1+1) + (1-R.thisq)*(utility(nd+1, ng+1)));
                
else
    
    QD = 0;
    
end

Qsa = [QG; QB; QD];