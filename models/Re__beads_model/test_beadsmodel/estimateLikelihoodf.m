function [logLikelihood, pickTrial, dQvec, ddec, aQvec, choice] = estimateLikelihoodf(alpha,Cw,q,Cs, sequence, aqvec_switch)

% this will be used for stopping at optimal draws (position)
findPick        = 1;

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
            numGreen    = numGreen + 1; % update numGreen  
        end
        
        % also increment draws 
        numDraws        = numDraws + 1;
        
        % compute action values for each draw (in current sequence)
        [v, d, Qvec]    = Val(q, numDraws, numGreen, alpha, sequenceL, Cw, Cs);
        
        % append actios values of the current sequence to dQvec
        dQvec(j,:)      = Qvec;
        
        % append choice probabilities 
        ddec(j,:)       = d;
        
        aQvec{i}(j,:)   = Qvec;
        
        % determine optimal stopping position 
        if findPick == 1 & j < nchoices & (Qvec(1) > Qvec(3) | Qvec(2) > Qvec(3))
            pickTrial(i) = j;
            break
            
        elseif findPick == 1 & j == nchoices
            pickTrial(i) = j;
        end
        
        % 
        if j == sequenceL
            d = [d; 0];
        end   
        
    end % end of draws loop
    
    % accumulate chosen urns
    [biggest_value choice(i)] = max(Qvec);

end % end of trials loop

end