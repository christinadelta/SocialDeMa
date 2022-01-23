function utility_t = stateUtilityBeads(utility, i, draw, maxDraws, numGreen, R)

utility_t = zeros(maxDraws, maxDraws + 1);

futureDraws = i - draw;

ndf = i;

for greenDraw = 0 : futureDraws
    
    ngf = numGreen + greenDraw;
    
    Qsa = actionValueBeads(utility, R, ndf, ngf, i, maxDraws);
    
    utility_t(ndf, ngf+1) = max(Qsa);
    
end % end of for loop
end