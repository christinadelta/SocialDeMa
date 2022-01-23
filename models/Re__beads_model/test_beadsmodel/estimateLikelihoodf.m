function [logLikelihood, pickTrial, dQvec, ddec, aQvec choice] = estimateLikelihoodf(alpha,Cw,q,Cs, sequence, aqvec_switch)

% ?
findPick1       = 1;

sequence        = seq_mat;              % rename sequence 
sequenceL       = size(sequence, 2);    % length of sequence 
ntrials         = size(sequence, 1);    % number of trials/sequences

logLikelihood   = 0;                    % initialise likelihood log to zero 

% loop over trials 
for i = 1:ntrials
    
    % chooce sequence
    thisSequence    = sequence(i,:);
    
    % number of choices for this sequence
    nchoices        = sequenceL;
    
    % initialise number of draws and number green beads (zero)
    numDraws        = 0;
    numGreen        = 0;
    
    dQvec           = []; % values for each action (vector)
    ddec            = []; % corresponding probabilities generated with softmax and alpha 
    
    % loop over draws for that sequence 
    for j = 1:nchoices
        
        % green or blue bead?
        if thisSequence(j) == 1
            numGreen = numGreen + 1; % update numGreen  
        end
        
        % increment draws 
        numDraws    = numDraws + 1;
        
        % compute action values for each draw (in current sequence)
        [v, d, Qvec] = Val(q, numDraws, numGreen, alpha, sequenceL, Cw, Cs)
        
    end % end of draws loop

end % end of trials loop



end