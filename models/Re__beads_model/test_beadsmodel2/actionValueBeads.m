function Qsa = actionValueBeads(utility, R, nd, ng, i, maxDraws)

pg = PG(R.q, nd, ng);

pb = 1 - pg;

QG = R.correct * pg + R.error * pg;
QB = R.correct * pb + R.error * pg;

if i < maxDraws
    
    QD = R.sample + pb*((1-R.q)*utility(nd+1, ng+1+1) +   (R.q)*(utility(nd+1, ng+1))) + ...
                    pg*(  (R.q)*utility(nd+1, ng+1+1) + (1-R.q)*(utility(nd+1, ng+1)));
                
else
    
    QD = 0;
                
end 

Qsa = [QG; QB; QD];

end