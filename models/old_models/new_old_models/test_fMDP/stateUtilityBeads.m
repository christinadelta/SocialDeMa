function utility_t = stateUtilityBeads(utility, drawi, draw, maxdraws, ng, Cc, Cw, Cs, thisq)

utility_t = zeros(maxdraws, maxdraws+1);

futureDraws = drawi - draw;

ndf = drawi;

for greenDraws = 0 : futureDraws
    
    ngf = ng + greenDraws;

    Qsa = actionValueBeads(utility, Cc, Cw, Cs, thisq, ndf, ngf, drawi, maxdraws);

    utility_t(ndf, ngf+1) = max(Qsa);        
    
end
return