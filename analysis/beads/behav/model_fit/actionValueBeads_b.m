function Qsa = actionValueBeads_b(utility, R, nd, ng, drawi, maxDraws)

pg = PG_b(R.thisq, nd, ng);

pb = 1 - pg;

% compute action values for each urn either by using the cost for being
% incorrect + reward for being correct, or by using just the difference
% between the two
% if R.model == 4
%     QG = R.diff * pb;
%     QB = R.diff * pg;
% else
    QG = R.costreward * pg + R.costloss * pb;
    QB = R.costreward * pb + R.costloss * pg;
% end

if drawi < maxDraws

    QD = R.sample + pb*((1-R.thisq)*utility(nd+1, ng+1+1) +   (R.thisq)*(utility(nd+1, ng+1))) + ...
                    pg*(  (R.thisq)*utility(nd+1, ng+1+1) + (1-R.thisq)*(utility(nd+1, ng+1)));
                
else
    
    QD = 0;
    
end

Qsa = [QG; QB; QD];

end