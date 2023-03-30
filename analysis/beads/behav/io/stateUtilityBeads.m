function utility_t = stateUtilityBeads(utility, drawi, draw, maxDraws, ng, R)

utility_t = zeros(maxDraws, maxDraws+1);

futureDraws = drawi - draw;

ndf = drawi;

for greenDraws = 0 : futureDraws
    
    ngf = ng + greenDraws;

    Qsa = actionValueBeads(utility, R, ndf, ngf, drawi, maxDraws);

    utility_t(ndf, ngf+1) = max(Qsa);        
    
end

return