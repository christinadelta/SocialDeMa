%%% this test script will be used to fitt Bruno's model to participant (or
%%% now simulated data) 

param = -0.25;
for cond = 1:conds 
    param = -0.25;
    % first deal with sequences 
    % which sequences to use? which condition are we in?
    if cond == 1
        seq_mat = sequence_e;
        condraws = drawse;
    else
        seq_mat = sequence_d;
        condraws = drawsd;
    end
    
    % get this cond urns
    condurn     = urns(:,cond);
    
    for len = 1:length(condurn)
        
        % deal with sequences first
        % 1. if green urn switch indecies in the sequence
        if condurn(len) == 0
            seq_ones = find( seq_mat(len,:) == 1);
            seq_twos = find( seq_mat(len,:) == 2);
            seq_mat(len,seq_ones) = 2;
            seq_mat(len,seq_twos) = 1;
        end

    end
    
    % this is for the model
    seq_mat(find(seq_mat==2)) = 0;
    
    % add thi condition's probability to the R struct
    R.q = probs(1,cond);
    R.cond = cond;
    
    options         = optimset('MaxFunEvals', 5000, 'TolFun', 0.001);
    
    llaMin          = Inf;
    startParam      = param;
    
    [mparams, lla] = fminsearch(@(param) fitmdp_beadsb(param, R, seq_mat, choiceVec, condraws),startParam, options);
    
    if lla < llaMin
        llaMin = lla;
        minParams = mparams;
    end
    
    fprintf('ll %.3f\n', lla);
    fprintf('min params %.3f\n', minParams);

    
    [ll, Qsad, cprob] = fitmdp_beads(minParams, R, seq_mat, choiceVec, condraws);
 
end % end of condition loop