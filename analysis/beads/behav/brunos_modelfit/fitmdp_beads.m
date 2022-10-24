function [ll, Qsad, cprob] = fitmdp_beads(param, R, thiscond_seqmat, thiscond_choiceVec)

maxDraws    = size(thiscond_seqmat,2); 
k           = 3;
R.sample    = param(1);
%beta        = param(2); % 
beta        = R.initbeta; %

% get number of trials/sequences 
ntrials     = max(size(thiscond_seqmat));
Qsad        = nan(ntrials,maxDraws,k);

% init log-likelihood
ll          = 0;

for trial = 1:ntrials
    
    % extract choice vector for this trial
    this_choices    = thiscond_choiceVec{1,trial};
    
    % what was the number of draws for this trial/sequence
    tdraws          = size(this_choices,1);

    % ok, now based on the number of draws for this trial, extract
    % the beads colours
    trialdraws      = thiscond_seqmat(trial,1:tdraws); 
    ndraws          = length(trialdraws); % this trial length (this is actually the same as number of draws)

    for draw = 1:ndraws  
        % run backwardutility
        Qsad(trial,draw,1:3)    = backWardUtility_b(trialdraws,draw,maxDraws,R)';
        vVec                    = Qsad(trial,draw,1:3);
        cprob(trial,draw,:)     = exp(beta*vVec)./sum(exp(beta*vVec));
     
    end % end of draws loop
    
    % not sure what the next line of code does
    % seqChoice = (choiceVec{cond,trial}(end) == 2)

    % so.... seqChoice should be either 2 or 1. In Bruno's code if
    % green then seqChoice = 1, if blue then seqChoice = 2
    % I think this is how it works

    if this_choices(end,1) == 1
        % then it's a blue trial (1)
        seqChoice = 1;
    elseif this_choices(end,2) == 1
        % then it's a green trial (2)
        seqChoice = 2;
    end

    if ndraws - 1 > 0 & ndraws < maxDraws
        ll = ll - sum(log(squeeze(cprob(trial, ndraws-1, k)))) - log(squeeze(cprob(trial, ndraws, seqChoice)));
    elseif ndraws < maxDraws
        ll = ll - log(squeeze(cprob(trial, ndraws, seqChoice)));
    end % end of if

    % llall(trial,cond) = ll;
    
end % end of trials loop

return