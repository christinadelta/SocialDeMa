function utCont = computeContinue(utility, postProb, x, ti)

postProb.nu = ti-1;

utCont = zeros(length(x), 10);

% pspx = zeros(length(x), length(x));

expData.n   = 1;
expData.sig = 0;

for xi = 1 : length(x)
    
    expData.mu  = x(xi);
    
    postProb = normInvChi(postProb, expData);
    spx = posteriorPredictive(x, postProb);
    spx = (spx/sum(spx));
    
    %     pspx(:, xi) = spx;
    
    % I am not sure if this is correct. should it be 100x1 or 100x10?
    utCont(xi) = spx'*utility;
    % utCont(xi,:) = spx'*utility;
    
end

reurn