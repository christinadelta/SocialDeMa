function [ll, Qsad, cprob] = fitmdp_beads(param, R, thiscond_seqmat, cond_choices,fitm)

maxDraws    = size(thiscond_seqmat,2); 
k           = 3;

% which parameters are free?
if R.model == 1
    R.sample        = R.initsample; % fixed
    beta            = param;        % free 
    R.costloss      = R.error;      % fixed 
    R.costreward    = R.correct;    % fixed 
elseif R.model == 2
    R.sample        = param(1);
    R.costloss      = R.error;
    R.costreward    = R.correct;
    beta            = param(2);
end 


% get number of trials/sequences 
ntrials         = max(size(thiscond_seqmat));
Qsad            = nan(ntrials,maxDraws,k);

% init log-likelihood
ll              = 0;

for trial = 1:ntrials
    
    % if this fuction runs to fitt the model and compute -ll then use
    % participant data
    if fitm == 1
        this_choices    = cond_choices{1,trial}; % extract choice vector for this trial
        tdraws          = size(this_choices,1);  % what was the number of draws for this trial/sequence

    else % use the maximum number of draws
        tdraws = maxDraws;

    end

    % ok, now based on the number of draws for this trial, extract
    % the beads colours from sequence mat
    trialdraws      = thiscond_seqmat(trial,1:tdraws); 
    ndraws          = length(trialdraws); % this trial length (this is actually the same as number of draws)

    for draw = 1:ndraws  
        % run backwardutility
        Qsad(trial,draw,1:3)    = backWardUtility_b(trialdraws,draw,maxDraws,R)';
        vVec                    = squeeze(Qsad(trial,draw,:))';
       
        cprobVals               = exp(beta*vVec)./sum(exp(beta*vVec)); % softmax to convert action values to probabiltiies of these actions
        cvec                    = cprobVals';
        cprob(trial,draw,:)     = cprobVals; 
        
        % if model fitting:
        if fitm == 1
            ll                  = ll - log(this_choices(draw,:)*cvec);
        end

%         if isnan(ll)
% 
%             fprintf('')
%         end
% 
%         if isnan(any(cvec))
%             fprintf('')
%         end

     
    end % end of draws loop
    
end % end of trials loop

% for i = 1:size(cvec,2)
%     fprintf(' choice probabilities %3.3f\n',cvec(i))
% end

fprintf(' beta %3.3f\n',beta)
fprintf(' ll %3.3f\n',ll)


return