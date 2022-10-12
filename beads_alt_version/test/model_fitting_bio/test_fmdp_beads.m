
%%% this test script will be used to fitt Bruno's model to participant (or
%%% now simulated data) 

beta = 0.13; % not sure if will use that (Nick uses ailpha value of 1)
alpha = R.alpha;
maxDraws = 10; 
k = 3;

for cond = 1:conds 
    
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
    
    % get number of trials/sequences 
    ntrials     = max(size(seq_mat));
    Qsad        = nan(ntrials,maxDraws,k);
    
    % add thi condition's probability to the R struct
    R.q = probs(1,cond);
    
    % init log-likelihood
    ll          = 0;
    
    for trial = 1:ntrials
        
        % what was the number of draws for this trial/sequence
        tdraws = condraws(trial);
        
        % ok, now based on the number of draws for this trial, extract
        % the beads colours
        trialdraws = seq_mat(trial,1:tdraws); % e.g. 2 draws, only both green (because it's a green urn trial)
        ndraws = length(trialdraws); % this trial length (this is actually the same as number of draws)
        
        for draw = 1:ndraws
            
            % run backwardutility
            Qsad(trial,draw,1:3) = backWardUtility_b(trialdraws,draw,maxDraws,R)'
            vVec = Qsad(trial,draw,1:3);
            cprob(trial,draw,:) = exp(beta*vVec)./sum(exp(beta*vVec))
     
        end % end of draws loop
        
        % not sure what the next line of code does
        % seqChoice = (choiceVec{cond,trial}(end) == 2)
        
        % so.... seqChoice should be either 2 or 1. In Bruno's code if
        % green then seqChoice = 1, if blue then seqChoice = 2
        % I think this is how it works
        
        if choiceVec{cond,trial}(end,1) == 1
            % then it's a blue trial (1)
            seqChoice = 1;
        elseif choiceVec{cond,trial}(end,2) == 1
            % then it's a green trial (2)
            seqChoice = 2;
        end
        
        if ndraws - 1 > 0 & ndraws < maxDraws
            ll = ll - sum(log(squeeze(cprob(trial, ndraws-1, k)))) - log(squeeze(cprob(trial, ndraws, seqChoice)))
        elseif ndraws < maxDraws
            ll = ll - log(squeeze(cprob(trial, ndraws, seqChoice)))
        end % end of if
       
    end % end of trials loop
  
end % end of condition loop