

    %%% not used anymore 
    
    options         = optimset('Display','iter','MaxFunEvals', 5000, 'TolFun', 0.001,'PlotFcns', @optimplotfval);
    
    llaMin          = Inf;
    startParam      = param;
    
    [mparams, lla] = fminsearch(@(param) fitmdp_beadsb(param, R, seq_mat, choiceVec, condraws),startParam, options);
    
    if lla < llaMin
        llaMin = lla;
        minParams = mparams;
    end
    
    fprintf('ll %.3f\n', lla);
    fprintf('optimal cost-sample: %.3f\n optimal beta: %.3f\n', minParams);

    
    [ll, Qsad, cprob] = fitmdp_beads(minParams, R, seq_mat, choiceVec, condraws);
 
