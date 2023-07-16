function [simR,simoutput] = simBeadsData(simvars, simR)

% created 31/03/23 as part of the BEADS task in OPTIMAL STOPPING PROBLEMS

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

% ------------

% unpack sim structure
ntrials         = simvars.ntrials;
cond_trls       = simvars.contrials;
maxDraws        = simvars.maxDraws;
qvals           = simvars.qvals;
conditions      = simvars.conditions;
cond            = simR.cond;
simR.k          = 3;

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
    % tmp_sequences(trl,:)            = (rand(1, maxDraws) > thisq) + 1; 

    clear s s_random 
end % end of trials loop

sim_sequnces                        = tmp_sequences;
drawSequence                        = tmp_sequences;

simR.thisq                          = thisq;
% simR.sample                         = simR.Cs; 

% deal with simulated urns, we need that to get simulated agent's
% choices 
for u = 1:length(urntype)
    if urntype(u) == 0 % if green urn switch index coding
        seq_ones = find(drawSequence(u,:) == 1);
        seq_twos = find(drawSequence(u,:) == 2);
        drawSequence(u,seq_ones) = 2;
        drawSequence(u,seq_twos) = 1;
    end 
end

% recode 2s to 0s for backward induction 
drawSequence(find(drawSequence==2))      = 0;

% get action values and simulate responses using softmax 
[~, Qsat, cprob] = mdp_beads(simR, drawSequence);

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

end % end of trials loop

simoutput.performance   = mean(correctResp);
simoutput.avsamples     = mean(cprob_samples);
simoutput.simsequences  = drawSequence;
simoutput.simurns       = simurns;
simoutput.simdraws      = cprob_samples;
simoutput.simchoice     = correctResp;
simoutput.simchoicevec  = simchoiceVectors;

end % end of function