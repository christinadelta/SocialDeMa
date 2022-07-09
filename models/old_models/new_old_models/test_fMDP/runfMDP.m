%% run bruno's fMDP model with pilot data

beta        = 0.13; % not sure the need for this
maxdraws    = 10; % number of max draws
ntrials     = max(size(cond_sequence)); % number of sequences 
Qsad        = NaN*ones(ntrials, maxdraws, 3);
ll          = 0; %init ll
sample      = Cs;
Cc          = 10;

% extract info
thisq       = info.p;
thisurns    = info.urntypes;
numdraws    = info.numdraws;

for trial = 1:ntrials
    
    trialdraws      = numdraws(trial);
    trialchoices    = cond_choices{1,trial};
    sequence        = cond_sequence{1,trial};
    
    % loop over draws 
    for draw = 1:maxdraws
        Qsad(trial,draw, 1:3)   = backWardUtility(sequence, draw, maxdraws, Cc, Cw, Cs, thisq);
        
        vVec                    = Qsad(trial,draw, 1:3);
        cprob(trial,draw,:)     = exp(beta*vVec)./sum(exp(beta*vVec)); % this is d in Nick's code
        
    end % end of draw loop
    
    
end % end of trial loop