function [logLikelihood, pickTrial, dQvec, ddec, aQvec, choice] = estimateLikelihoodf(alpha,Cd,q,Cs, sequence, aqvec_switch)

% this will be used for stopping at optimal draws (position)
findPick            = 1;
lseq                = size(sequence, 2);    % length of sequence 
ntrials             = size(sequence, 1);    % number of trials/sequences

logLikelihood       = 0;                    % initialise likelihood log to zero 

% rename sequence
this_sequence        = sequence;

% number of choices for this sequence
nchoices            = lseq;

% initialise number of draws and number green beads (zero)
numDraws            = 0;
numGreen            = 0;

dQvec               = []; % values for each action (vector)
ddec                = []; % corresponding probabilities generated with softmax and alpha 

% loop over draws for that sequence 
for j = 1:nchoices

    % green or blue bead? or majority beads colour in the urn? -- I think it is
    % the majority colour 
    if this_sequence(j) == 1
        numGreen    = numGreen + 1; % update numGreen  
    end

    % also increment draws 
    numDraws        = numDraws + 1;

    % compute action values for each new draw (in current sequence)
    % until action value for one of the two urns exceeds action value
    % for drawing again. 
    [v, d, Qvec]    = Val(q, numDraws, numGreen, alpha, lseq, Cd, Cs);

    % append action values of the current sequence to dQvec
    dQvec(j,:)      = Qvec;

    % append choice probabilities 
    ddec(j,:)       = d;

    aQvec(j,:)      = Qvec;

    % determine optimal stopping position 
    if findPick == 1 & j < nchoices & (Qvec(1) > Qvec(3) | Qvec(2) > Qvec(3))
        pickTrial   = j; % number of draws
        break

    elseif findPick == 1 & j == nchoices
        pickTrial = j;
    end

    % 
    if j == lseq
        d = [d; 0];
    end   

end % end of draws loop

% accumulate chosen urns
[biggest_value choice] = max(Qvec);

end