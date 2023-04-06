function simoutput = simBeadsData(simvars, simR)

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
ntrials     = simvars.ntrials;
cond_trls   = ntrials/2;
maxDraws    = simvars.maxDraws;
qvals       = simvars.qvals;
conditions  = simvars.conditions;
k           = 3;

% simlute vector with urns 
temp            = [ones(cond_trls,1); zeros(cond_trls,1)]';
temp_rand       = randperm(length(temp));
temp            = temp(temp_rand);
simurns(1:cond_trls,1) = temp(1:cond_trls);
simurns(1:cond_trls,2) = temp(cond_trls+1:end);

% loop over conditions 
for cond = 1:conditions 
    
    % first simulate sequences 
    % extract probability value and urns for this condition
    thisq       = qvals(cond);
    urntype     = simurns(:,cond);

    % simulate sequences
    for trl = 1:cond_trls

        s                               = [ones(1,thisq*maxDraws), ones(1,maxDraws-thisq*maxDraws)*2]; % create tmp sequence
        s_random                        = randperm(length(s));
        tmp_sequences(trl,:)            = s(s_random); 

    clear s s_random 
    end % end of trials loop

    sim_sequnces{1,cond}               = tmp_sequences;
    drawSequence                       = tmp_sequences;

    % what is the probability of this condition? 
    simR.thisq                          = thisq;
    simR.sample                         = simR.initsample; 
    simR.beta                           = simR.initbeta;
    
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
        
    % run backward induction
    [r, Qsat]                           = backWardInduction(drawSequence, simR);

    % get responses and number of draws 
    for dri = 1:cond_trls

        choiceTrial             = find(squeeze(Qsat(dri, :, 3)) - max(squeeze(Qsat(dri, :, 1:2))') < 0);
        simDraws(dri,cond)      = choiceTrial(1);
        [ma ma_i]               = max(squeeze(Qsat(dri,simDraws(dri,cond),:))); % which of the two urn was chosen? [based on AQ values]

        % get simulated agent's choices 
        if (ma_i == 1 & urntype(dri) == 1) | (ma_i == 2 & urntype(dri) == 0)
            simChoice(dri,cond) = 1;
        else
            simChoice(dri,cond) = 0;
        end

        % simulate choice vectors 
        draws_choice                    = simDraws(dri,cond); 
        t                               = nan(draws_choice, k); % init empty matrix
        
        for d = 1:draws_choice
            
            if d ~= draws_choice % if this is not the last draw add 0's to b and g columns and 1 to draw column
                t(d,1:2)                = 0; % index zero for b and g columns
                t(d,3)                  = 1; % index one for draw column
            else
                if urntype(dri) == 1 & simChoice(dri,cond) == 1 % this is a blue urn and sub was correct
                    t(d,2:3)            = 0; % index zero for g and draw columns
                    t(d,1)              = 1; % index one for b column (sub ressed blue)
                elseif urntype(dri) == 1 & simChoice(dri,cond) == 0 % this is a blue urn and sub was incorrect or did not respond
                    t(d,1)              = 0; % index zero for b 
                    t(d,2)              = 1; % index one for g column
                    t(d,3)              = 0; % index zero for draw 
                elseif urntype(dri) == 0 & simChoice(dri,cond) == 1 % this is a green urn and sub was correct
                    t(d,1)              = 0; % index zero for b 
                    t(d,2)              = 1; % index one for g column
                    t(d,3)              = 0; % index zero for s 
                elseif urntype(dri) == 0 & simChoice(dri,cond) == 0 % this is a green urn and sub was incorrect or did not respond
                    t(d,2:3)            = 0; % index zero for g and draw columns
                    t(d,1)              = 1; % index one for b column   
                end
            end
        end % end of d draws loop

        % store choice vector
        simchoiceVectors{1,cond}{1,dri} = t;
        
    end % end of dri trials loop


    clear tmp_sequences
end % end of conditions loop


% compute mean performance of simulated choices
simoutput.performance   = mean(simChoice);
simoutput.avsamples     = mean(simDraws);

% sore all simulated output in struct
simoutput.simsequences  = sim_sequnces;
simoutput.simurns       = simurns;
simoutput.simdraws      = simDraws;
simoutput.simchoice     = simChoice;
simoutput.simchoicevec  = simchoiceVectors;

end % end of function