function [sim_sequnces, sim_choiceVectors] = simBeadsData(simvars)

% created 31/03/23 as part of the BEADS task in OPTIMAL STOPPING PROBLEMS

% Input: struct that contains all the variables needed to simulate sequences of bead draws
%       - number of trials/sequences
%       - maximu draws
%       - propaibility values (qvals)
%       - conditions 

% Outputs: 
%       - cell with 2 26-by-10 sequences with bead draws 
%       - cell with simulated choice vectors (10-by3)

% ------------

% unpack sim structure
ntrials     = simvars.ntrials;
cond_trls   = ntrials/2;
maxDraws    = simvars.maxDraws;
qvals       = simvars.qvals;
conditions  = simvars.conditions;
k           = 3;

% loop over conditions 
for cond = 1:conditions 

    % extract probability value for this condition
    thisq   = qvals(cond);

    for trl = 1:ntrials/2

        s                               = [ones(1,thisq*maxDraws), ones(1,maxDraws-thisq*maxDraws)*2]; % create tmp sequence
        s_random                        = randperm(length(s));
        tmp_sequences(trl,:)            = s(s_random); 

        % create a vector with choices for this sequence (should be up to
        % 10th draw)
        draw                            = maxDraws; % number of draws
        choicevec{trl}                  = zeros(draw, k); % set the matrix 
        choicevec{trl}((1:(end-1)),3)   = 1; % choices for drawing again up to last bead (ones in column 3)

        if mean(tmp_sequences(trl, 1:draw)) > 1.5
            choicevec{trl}(end, 2)      = 1;
        else
            choicevec{trl}(end, 1)      = 1;
        end

    clear s s_random 
    end % end of trials loop

    % recode 2s to 0s for backward induction 
     tmp_sequences(find(tmp_sequences==2))       = 0;

     sim_sequnces{1,cond}       = tmp_sequences;
     sim_choiceVectors{1,cond}  = choicevec;

     clear tmp_sequences choicevec
end % end of conditions loop


end % end of function