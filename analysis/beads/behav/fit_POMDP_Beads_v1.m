function modeloutput = fit_POMDP_Beads_v1(R, thiscond_model,best_params)


% comment will go here...
% and here..
 %% 

maxDraws    = size(thiscond_model.simsequences,2); 
k           = 3;
urntype     = thiscond_model.simurns;

% Initialize log-likelihood
ll              = 0;

% which parameters are free?
% Transform beta to be within [min_beta, max_beta]
beta            = best_params(2);
R.sample        = best_params(1);
% R.sample        = R.costsample;

R.costloss      = R.error;      % fixed 
R.costreward    = R.correct;    % fixed 


% get number of trials/sequences 
ntrials         = R.contrials;
Qsad            = nan(ntrials,maxDraws,k);
sequence_matrix = thiscond_model.simsequences;

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

% run backward induction
for trial = 1:ntrials
   
    this_choices    = thiscond_model.simchoicevec{1,trial}; % extract choice vector for this trial
    tdraws          = size(this_choices,1);                 % what was the number of draws for this trial/sequence
    
    % ok, now based on the number of draws for this trial, extract
    % the beads colours from sequence mat
    trialdraws      = sequence_matrix(trial,1:tdraws); 
    ndraws          = length(trialdraws);                   % this trial length (this is actually the same as number of draws)

    for draw = 1:ndraws  
        % run backwardutility
        Qsad(trial,draw,1:3)    = backWardUtility(trialdraws,draw,maxDraws,R)';
        vVec                    = squeeze(Qsad(trial,draw,:))';
       
        cprobVals               = exp(beta*vVec)./sum(exp(beta*vVec)); % softmax to convert action values to probabiltiies of these actions
        cvec                    = cprobVals';
        cprob(trial,draw,:)     = cprobVals; 
        
    end % end of draws loop   
end % end of trials loop

% use the choice probabilities to compute choices 
N                               = 1000;
[cprob_samples,model_urnchoice] = modelSamples(cprob,N);

% store model output structure
modeloutput.samples             = cprob_samples;
modeloutput.actionVals          = Qsad;
modeloutput.avsamples           = mean(cprob_samples);

end % end of main function

%% run backward utility 

function Qsa = backWardUtility(drawSequence, draw, maxDraws, R)

utility = zeros(maxDraws, maxDraws+1);

ng = sum(drawSequence(1:draw));

for drawi = maxDraws : -1 : (draw + 1)
        
    [utility] = stateUtilityBeads(utility, drawi, draw, maxDraws, ng, R);
    
end
    
Qsa = actionValuesBeads(utility, R, draw, ng, draw, maxDraws);
    
end % end of backward utility function

%% run state utility and generate action values

function utility_t = stateUtilityBeads(utility, drawi, draw, maxDraws, ng, R)

utility_t = zeros(maxDraws, maxDraws+1);

futureDraws = drawi - draw;

ndf = drawi;

for greenDraws = 0 : futureDraws
    
    ngf = ng + greenDraws;

    Qsa = actionValuesBeads(utility, R, ndf, ngf, drawi, maxDraws);

    utility_t(ndf, ngf+1) = max(Qsa);        
    
end

end % end of state utility function

function Qsa = actionValuesBeads(utility, R, nd, ng, drawi, maxDraws)

pg = PG(R.thisq, nd, ng);

pb = 1 - pg;

QG = R.correct*pg + R.error*pb;
QB = R.correct*pb + R.error*pg;

if drawi < maxDraws

    QD = R.sample + pb*((1-R.thisq)*utility(nd+1, ng+1+1) + (R.thisq)*(utility(nd+1, ng+1))) + ...
                    pg*((R.thisq)*utility(nd+1, ng+1+1) + (1-R.thisq)*(utility(nd+1, ng+1)));
                
else
    
    QD = 0;
    
end

Qsa = [QG; QB; QD];

end % end of action values function

%% compute probabilites function 
function p = PG(q, nd, ng)

p = 1/(1 + (q/(1-q))^(nd-2*ng));

end % end of probabilities function

%% compute model samples 

function [cprob_samples,model_urnchoice] = modelSamples(cprob,N)

% compute model sampling rate using cprob instead
% created 24/03/2023

% ---------------------------

len     = 10;                   % true sequence length
trials  = length(cprob(:,1,1)); % number of sequences

% if len 

% loop over sequences 
for t = 1:trials

    choiceProbs                 = squeeze(cprob(t,:,:)); % extract this trial choice probabilities 

    % add a stopping point at the competing urn at last draw (just to
    % ensure that model stops drawing)
    [val urnchoice]             = max(choiceProbs(end,1:2));
    choiceProbs(end,urnchoice)  = Inf;
    model_urnchoice(t,1)        = urnchoice;

    % sometimes the model draws les than the actual sequence length; find rows with sum zero...
    for l = 1:size(choiceProbs,1)

        if choiceProbs(l,:) == sum(0)
            choiceProbs(l,:) = nan;
        end


    end % end of cprob length loop
    % ...and remove them
    
    choiceProbs(any(isnan(choiceProbs), 2), :) = [];

    % in cases that the model uses less beads than the total sequence (e.g., 6 draws instaed of 10), the
    % sequence vector we get is 6 rows long; I fill the rest of column 3
    % with -inf and with 1s the column of the competing/chosen urn 
    % THIS PART IS NOT REALLY NEEDED FOR NOW
    if size(choiceProbs,1) < len

        if urnchoice == 1 % if model chose the blue urn
            tmp_diff                        = [ones(len,1) zeros(len,1) zeros(len,1)]; % [1 0 0]
        elseif urnchoice == 2 % if model chose the green urn
            tmp_diff                        = [zeros(len,1) ones(len,1) zeros(len,1)]; % [0 1 0]
        end

        tmp_diff(1:size(choiceProbs,1),:)   = choiceProbs;
        choiceProbs                         = tmp_diff;
    end 

    % ok now its time to compute model's sampling rates based on the 
    for i = 1:N

        isample             = rand(1,size(choiceProbs,1))';

        % now that we have this information let's see what the model is choosing based on the "softmax" choice probabilities 
        tmp_samples(i,1)    = find(choiceProbs(:,urnchoice) > isample(:,1),1, "first");
        clear isample
    end % end of sampling loop

    cprob_samples(t,1)      = round(mean(tmp_samples)); % round to nearest decimal

    clear choiceProbs tmp_samples

end % end of trials loop

end % end of model samples function
