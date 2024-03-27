function simoutput = runIdealObserver(simvars)

% main function 
% created in March 2024 as part of the BEADS task in OPTIMAL STOPPING PROBLEMS

% Input: struct that contains all the variables needed to simulate sequences of bead draws
%       - number of trials/sequences
%       - maximu draws
%       - propaibility values (qvals)
%       - conditions 

% Outputs: 
%       - structure that contains:
%                   - simulated choices, performance number of draws
%                   - simulated choice vectors and sequences
%                   - simulated urns

%% unpack simvar structure

ntrials         = simvars.ntrials;
cond_trls       = simvars.contrials;
maxDraws        = simvars.maxDraws;
qvals           = simvars.qvals;
conditions      = simvars.conditions;
cond            = simvars.cond;
simvars.k       = 3;

%% simulate some sequences to work with

% simulate condition urns
temp            = [ones(cond_trls/2,1); zeros(cond_trls/2,1)]';
temp_rand       = randperm(length(temp));
temp            = temp(temp_rand);
simurns         = temp;

% first simulate sequences 
% extract probability value and urns for this condition
thisq           = qvals(cond);
urntype         = simurns;

% simulate sequences
for trl = 1:cond_trls

    s                               = [ones(1,thisq*maxDraws), ones(1,maxDraws-thisq*maxDraws)*2]; % create tmp sequence
    s_random                        = randperm(length(s));
    tmp_sequences(trl,:)            = s(s_random);  

    clear s s_random 
end % end of trials loop

sim_sequnces                        = tmp_sequences;
drawSequence                        = tmp_sequences;
simSequence                         = tmp_sequences; % for saving (this one is coded as 1s and 2s)
simvars.thisq                       = thisq;

% deal with simulated urns, we need that to get simulated agent's
% choices 
for u = 1:length(urntype)
    if urntype(u) == 0 % if green urn switch index coding
        seq_ones = find(drawSequence(u,:) == 1);
        seq_twos = find(drawSequence(u,:) == 2);
        drawSequence(u,seq_ones)    = 2;
        drawSequence(u,seq_twos)    = 1;
    end 
end

% recode 2s to 0s for backward induction 
drawSequence(find(drawSequence==2)) = 0;

% get action values and simulate responses using softmax 
[~, Qsat, cprob] = mdp_beads(simvars, drawSequence);

% use the choice probabilities to compute choices 
N                               = 1000;
[cprob_samples,model_urnchoice] = modelSamples(cprob,N);

% get performance and choice vectors based on num samples and urn_choices 
for dri = 1:size(model_urnchoice,1)          
    if (model_urnchoice(dri) == 1 & urntype(dri) == 1) | (model_urnchoice(dri) == 2 & urntype(dri) == 0)
        correctResp(dri) = 1;
    else
        correctResp(dri) = 0;
    end % end of condition
    
    % simulate choice vectors
    draws_choice = cprob_samples(dri);

    if draws_choice == 1
        draws_choice    = draws_choice + 1; % to allow the model to use one more sample as in the task participants had to at least sample once after the first bead 
    end
    t            = nan(draws_choice, 3); % init empty matrix

    for d = 1:draws_choice
        
        if d ~= draws_choice % if this is not the last draw add 0's to b and g columns and 1 to draw column
            t(d,1:2)                = 0; % index zero for b and g columns
            t(d,3)                  = 1; % index one for draw column
        else
            if urntype(dri) == 1 & correctResp(dri) == 1 % this is a blue urn and sub was correct
                t(d,2:3)            = 0; % index zero for g and draw columns
                t(d,1)              = 1; % index one for b column (sub ressed blue)
            elseif urntype(dri) == 1 & correctResp(dri) == 0 % this is a blue urn and sub was incorrect or did not respond
                t(d,1)              = 0; % index zero for b 
                t(d,2)              = 1; % index one for g column
                t(d,3)              = 0; % index zero for draw 
            elseif urntype(dri) == 0 & correctResp(dri) == 1 % this is a green urn and sub was correct
                t(d,1)              = 0; % index zero for b 
                t(d,2)              = 1; % index one for g column
                t(d,3)              = 0; % index zero for s 
            elseif urntype(dri) == 0 & correctResp(dri) == 0 % this is a green urn and sub was incorrect or did not respond
                t(d,2:3)            = 0; % index zero for g and draw columns
                t(d,1)              = 1; % index one for b column   
            end
        end
    end % end of d draws loop

    % store choice vector
    simchoiceVectors{1,dri} = t;

    simoutput.performance   = mean(correctResp);
    simoutput.avsamples     = mean(cprob_samples);
    simoutput.simsequences  = simSequence;
    simoutput.simurns       = simurns;
    simoutput.simdraws      = cprob_samples;
    simoutput.simchoice     = correctResp;
    simoutput.simchoicevec  = simchoiceVectors;

end % end of trials loop

end % end of function 

%% mdb_beads function to run backward utility and softmax

function [ll, Qsat, cprob] = mdp_beads(simvars, drawSequence)

maxDraws            = size(drawSequence,2); 
k                   = simvars.k;
beta                = simvars.beta;
simvars.sample      = simvars.costsample;


% get number of trials/sequences 
ntrials         = max(size(drawSequence));
Qsat            = nan(ntrials,maxDraws,k);

% init log-likelihood
ll              = 0;

for trial = 1:ntrials

    Qsad            = zeros(maxDraws, k); % action values for this sequence will be stored here
    thisSequence    = drawSequence(trial,:);
    
    for draw = 1 : maxDraws
                   
        Qsad(draw, :) = backWardUtil(thisSequence, draw, maxDraws, simvars);
        vVec          = Qsad(draw,:);

        % compute choices using softmax 
        cprobVals               = exp(beta*vVec)./sum(exp(beta*vVec)); % softmax to convert action values to probabiltiies of these actions
        cvec                    = cprobVals';
        cprob(trial,draw,:)     = cprobVals; 
        
    end

    % store this_trial Q values in a cell
    Qsat(trial,:,:) = Qsad;

    clear Qsad 
end % end of trials loop

end % end of mdb_beads function

%% run backward utility 

function Qsa = backWardUtil(drawSequence, draw, maxDraws, simvars)

utility = zeros(maxDraws, maxDraws+1);

ng = sum(drawSequence(1:draw));

for drawi = maxDraws : -1 : (draw + 1)
        
    [utility] = stateUtilBeads(utility, drawi, draw, maxDraws, ng, simvars);
    
end
    
Qsa = compActionValues(utility, simvars, draw, ng, draw, maxDraws);
    
end % end of backward utility function

%% run state utility and generate action values

function utility_t = stateUtilBeads(utility, drawi, draw, maxDraws, ng, simvars)

utility_t = zeros(maxDraws, maxDraws+1);

futureDraws = drawi - draw;

ndf = drawi;

for greenDraws = 0 : futureDraws
    
    ngf = ng + greenDraws;

    Qsa = compActionValues(utility, simvars, ndf, ngf, drawi, maxDraws);

    utility_t(ndf, ngf+1) = max(Qsa);        
    
end

end % end of state utility function

function Qsa = compActionValues(utility, simvars, nd, ng, drawi, maxDraws)

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
