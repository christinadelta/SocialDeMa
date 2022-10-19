function [ll, Qsad, cprob] = fitmdp_beads(param, R, seq_mat, choiceVec, condraws)


% beta = 0.13; % not sure if will use that (Nick uses ailpha value of 1)

% alpha = R.alpha;
maxDraws    = 10; 
k           = 3;
R.sample    = param(1);
R.beta      = param(2);

% get number of trials/sequences 
ntrials     = max(size(seq_mat));
Qsad        = nan(ntrials,maxDraws,k);

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
        Qsad(trial,draw,1:3) = backWardUtility_b(trialdraws,draw,maxDraws,R)';
        vVec = Qsad(trial,draw,1:3);
        cprob(trial,draw,:) = exp(R.beta*vVec)./sum(exp(R.beta*vVec));
     
    end % end of draws loop
    
    % not sure what the next line of code does
    % seqChoice = (choiceVec{cond,trial}(end) == 2)

    % so.... seqChoice should be either 2 or 1. In Bruno's code if
    % green then seqChoice = 1, if blue then seqChoice = 2
    % I think this is how it works

    if choiceVec{R.cond,trial}(end,1) == 1
        % then it's a blue trial (1)
        seqChoice = 1;
    elseif choiceVec{R.cond,trial}(end,2) == 1
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