function ll = beads_lik_v1(param, simvars, sequence_matrix, condition_choices)

maxDraws    = size(sequence_matrix,2); 
k           = 3;
urntype     = simvars.urntype;

% Initialize log-likelihood
ll              = 0;

% parameters to be estimated
simvars.sample          = param(1);
% simvars.sample          = simvars.costsample;
beta                    = param(2);

% fixed parameters 
simvars.costloss      = simvars.error;
simvars.costreward    = simvars.correct;

% get number of trials/sequences 
ntrials         = max(size(sequence_matrix));
Qsad            = nan(ntrials,maxDraws,k);

% deal with simulated urns, we need that to get simulated agent's
% choices 
for u = 1:length(urntype)
    if urntype(1,u) == 0 % if green urn switch index coding
        seq_ones = find(sequence_matrix(u,:) == 1);
        seq_twos = find(sequence_matrix(u,:) == 2);
        sequence_matrix(u,seq_ones)    = 2;
        sequence_matrix(u,seq_twos)    = 1;
    end 
end

% recode 2s to 0s for backward induction 
sequence_matrix(find(sequence_matrix==2)) = 0;

for trial = 1:ntrials

    % extract this trial's sequence and agent's choices
    this_choices    = condition_choices{1,trial}; % extract choice vector for this trial
    tdraws          = size(this_choices,1);  % what was the number of draws for this trial/sequence

    % ok, now based on the number of draws for this trial, extract
    % the beads colours from sequence mat
    trialdraws      = sequence_matrix(trial,1:tdraws); 
    ndraws          = length(trialdraws); % this trial length (this is actually the same as number of draws)

    for draw = 1:ndraws  
        % run backwardutility
        Qsad(trial,draw,1:3)    = backWardUtility(trialdraws,draw,maxDraws,simvars)';
        vVec                    = squeeze(Qsad(trial,draw,:))';
       
        cprobVals               = exp(beta*vVec)./sum(exp(beta*vVec)); % softmax to convert action values to probabiltiies of these actions
        cvec                    = cprobVals';
        cprob(trial,draw,:)     = cprobVals; 
        ll                      = ll - log(this_choices(draw,:)*cvec);

    end % end of draws loop
end % end of trials loop
end % end of main function 

%% run backward utility 

function Qsa = backWardUtility(drawSequence, draw, maxDraws, simvars)

utility = zeros(maxDraws, maxDraws+1);

ng = sum(drawSequence(1:draw));

for drawi = maxDraws : -1 : (draw + 1)
        
    [utility] = stateUtilityBeads(utility, drawi, draw, maxDraws, ng, simvars);
    
end
    
Qsa = actionValuesBeads(utility, simvars, draw, ng, draw, maxDraws);
    
end % end of backward utility function

%% run state utility and generate action values

function utility_t = stateUtilityBeads(utility, drawi, draw, maxDraws, ng, simvars)

utility_t = zeros(maxDraws, maxDraws+1);

futureDraws = drawi - draw;

ndf = drawi;

for greenDraws = 0 : futureDraws
    
    ngf = ng + greenDraws;

    Qsa = actionValuesBeads(utility, simvars, ndf, ngf, drawi, maxDraws);

    utility_t(ndf, ngf+1) = max(Qsa);        
    
end

end % end of state utility function

function Qsa = actionValuesBeads(utility, simvars, nd, ng, drawi, maxDraws)

pg = PG(simvars.thisq, nd, ng);

pb = 1 - pg;

QG = simvars.correct*pg + simvars.error*pb;
QB = simvars.correct*pb + simvars.error*pg;

if drawi < maxDraws

    QD = simvars.sample + pb*((1-simvars.thisq)*utility(nd+1, ng+1+1) + (simvars.thisq)*(utility(nd+1, ng+1))) + ...
                    pg*((simvars.thisq)*utility(nd+1, ng+1+1) + (1-simvars.thisq)*(utility(nd+1, ng+1)));
                
else
    
    QD = 0;
    
end

Qsa = [QG; QB; QD];

end % end of action values function

%% compute probabilites function 
function p = PG(q, nd, ng)

p = 1/(1 + (q/(1-q))^(nd-2*ng));

end % end of probabilities function

