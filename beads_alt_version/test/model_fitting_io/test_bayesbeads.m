% script for fitting Nick's model on the data created in test_model.m

for cond = 1:2
    
    R.q = probs(cond);
    
    params = R.sample;
    fixedParams = [R.alpha; R.q; R.error; R.correct; cond];
    findPick = 1;
    condurns = urns(:,cond);
    
    % first deal with sequences 
    % which sequences to use? which condition are we in?
    if cond == 1
        seq_mat = sequence_e;
        condraws = drawse;
    else
        seq_mat = sequence_d;
        condraws = drawsd;
    end
    
    % % if green urn switch indecies 
    for s = 1:length(condurns)
        
        if condurns(s) == 0
            
            seq_ones = find( seq_mat(s,:) == 1);
            seq_twos = find( seq_mat(s,:) == 2);
            seq_mat(s,seq_ones) = 2;
            seq_mat(s,seq_twos) = 1;
            
        end
    end
    
    % this is for the model
    seq_mat(find(seq_mat==2)) = 0;
    
    Ntrials = size(seq_mat,1);
    maxDraws = size(seq_mat,2);
    
    options         = optimset('MaxFunEvals', 5000, 'TolFun', 0.001);
    
    llaMin          = Inf;
    
    startParam = params;
    
    [mparams, lla] = fminsearch(@(params) estimateLikelihood_n(params, seq_mat, choiceVec, fixedParams, findPick),startParam, options);
    
    if lla < llaMin
        llaMin = lla;
        minParams = mparams;
    end
    
    fprintf('ll %.3f\n', lla);
    fprintf('min params %.3f\n', minParams);

    % ftxt = sprintf('subjectParams_%d_%d.mat', subject,types);
    % save(ftxt, 'minParams', 'llaMin');

    % [ll, pickTrial, dQvec, ddec, aQvec] = estimateLikelihoodf(minParams, sequence, choiceVec, fixedParams, findPick);
    [ll, pickTrial, dQvec, ddec, aQvec] = estimateLikelihoodf_n(minParams, seq_mat, choiceVec, fixedParams, findPick);
    
    
end