function [ll, cprob, Qsad] = beadsmdp(params, thisub_choices, thisub_seq, info, R)

% extract alpha from fixed params struct (R)
alpha       = R.alpha;
Ntrials     = max(size(thisub_seq));
maxDraws    = 10;

Qsad        = NaN*ones(Ntrials, maxDraws, 3);

% init ll
ll          = 0;

for trial = 1 : Ntrials
    
    trialDraws                  = thisub_seq{1,trial};
    
    % convert twos in the sequence to zeros (otherwise calculation of majority 
    % bead colour will be wrong - see backwardUtility.m)
    seq_twos                    = find(trialDraws == 2);
    trialDraws(seq_twos)        = 0; 
     
    nDraws                      = length(trialDraws);
    
    R.sample                    = params(1);
    
    for draw = 1:nDraws
        
        Qsad(trial,draw, 1:3)   = backWardUtility(trialDraws,draw,maxDraws, R)';
        vVec                    = Qsad(trial,draw, 1:3);
        cprob(trial, draw, :)   = exp(alpha*vVec)./sum(exp(alpha*vVec));
        
    end % end of draws loop
    
    seqChoice = (thisub_choices{trial}(end) == 2) + 1; % ????
    
    if nDraws-1 > 0 & nDraws < maxDraws
        ll = ll - sum(log(squeeze(cprob(trial, nDraws-1, 3)))) - log(squeeze(cprob(trial, nDraws, seqChoice)));
%             ll = ll - (Qsad(trial, nDraws-1, 3) - max(Qsad(trial, nDraws-1, :))) - ...
%                       (Qsad(trial, nDraws, seqChoice) - Qsad(trial, nDraws, 3));
    elseif nDraws < maxDraws
%             ll = ll - (Qsad(trial, nDraws, seqChoice) - Qsad(trial, nDraws, 3));
        ll = ll - log(squeeze(cprob(trial, nDraws, seqChoice)));
    end  
    
end % end of trials loop 



return