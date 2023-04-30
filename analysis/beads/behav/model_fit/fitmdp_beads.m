function [ll, Qsad, cprob] = fitmdp_beads(param, R, thiscond_seqmat, cond_choices,fitm)

maxDraws    = size(thiscond_seqmat,2); 
k           = 3;

% which parameters are free?
if R.model == 1
    R.sample        = param;        % free
    beta            = R.beta;       % fixed 
    R.costloss      = R.error;      % fixed 
    R.costreward    = R.correct;    % fixed 
elseif R.model == 2
    R.costloss      = param(1);
    R.costreward    = param(2);
    R.sample        = R.Cs;
    beta            = R.beta;
elseif R.model == 3
    R.sample        = param(1);
    R.costloss      = param(2);
    R.costreward    = param(3);
    beta            = R.beta;
elseif R.model == 4
    R.diff          = param;
    R.sample        = R.Cs;
    beta            = R.beta;
elseif R.model == 5
    R.sample        = R.Cs;
    R.costloss      = R.error;
    R.costreward    = R.correct;
    beta            = param;
elseif R.model == 6
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
        
        % I'll test running softmax seperately for each of the three
        % options.. let's see if that works or helps with the beta
        % recovery!! 
%         cprob_blue              = exp(beta * vVec(1,1)) ./ (exp(beta * vVec(1,1))) + (exp(beta * vVec(1,2))) + (exp(beta * vVec(1,3)));
%         cprob_green             = exp(beta * vVec(1,2)) ./ (exp(beta * vVec(1,1))) + (exp(beta * vVec(1,2))) + (exp(beta * vVec(1,3)));
%         cprob_draw              = exp(beta * vVec(1,3)) ./ (exp(beta * vVec(1,1))) + (exp(beta * vVec(1,2))) + (exp(beta * vVec(1,3)));

        cprob_blue              = exp(beta * vVec(1,1)) ./ sum(exp(beta * vVec(1,:)));
        cprob_green             = exp(beta * vVec(1,2)) ./ sum(exp(beta * vVec(1,:)));
        cprob_draw              = exp(beta * vVec(1,3)) ./ sum(exp(beta * vVec(1,:)));
        cVec                    = [cprob_blue cprob_green cprob_draw]';
        
        % store all choice probabilities 
        cprob(trial,draw,1)     = cprob_blue;
        cprob(trial,draw,2)     = cprob_green;
        cprob(trial,draw,3)     = cprob_draw;
        
%         cprobVals               = exp(beta*vVec(1,1,:))./sum(exp(beta*vVec(1,1,:)));
%         %cprobVals               = exp(beta*vVec)./sum(exp(beta*vVec)); % softmax to convert action values to probabiltiies of these actions
%         cvec                    = squeeze(cprobVals);
%         cprob(trial,draw,:)     = cprobVals; 
        
        % if model fitting:
        if fitm == 1
            ll                  = ll - log(this_choices(draw,:)*cVec);
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