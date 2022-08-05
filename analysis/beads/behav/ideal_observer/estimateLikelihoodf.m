function [logLikelihood, pickTrial, dQvec, ddec, aQvec, choice] = estimateLikelihoodf(alpha,Cw, Cc, thisq,Cs, thiscond_seq, aqvec_switch)

% this will be used for stopping at optimal draws (position)
findPick            = 1; 
ntrials             = size(thiscond_seq, 2);    % number of trials/sequences
logLikelihood       = 0;                        % initialise likelihood log to zero  

% loop over condition sequences (26)
for trl = 1:ntrials
    
    % extract this trial's sequence
    thisequence = thiscond_seq{1,trl};
    
    % number of choices for this sequence
    lseq                = length(thisequence);
    nchoices            = lseq;
    
    % initialise number of draws and number green beads (zero)
    numDraws            = 0;
    numGreen            = 0;

    dQvec               = []; % values for each action (vector)
    ddec                = []; % corresponding probabilities generated with softmax and alpha 
    
    for draw = 1:nchoices 
        
        % if this_bead is of the majority colour: 
        if thisequence(draw) == 1
            numGreen    = numGreen + 1; % update numGreen  
        end

        % also increment draws 
        numDraws        = numDraws + 1;

        % compute action values for each new draw (in current sequence)
        % until action value for one of the two urns exceeds action value
        % for drawing again. 
        [v, d, Qvec]    = Val(thisq, numDraws, numGreen, alpha, lseq, Cw, Cc, Cs);
        
        % append action values of the current sequence to dQvec
        dQvec(draw, 1:length(Qvec))         = Qvec;

        % append choice probabilities 
        ddec(draw, 1:length(d))             = d;
        aQvec{trl}(draw, 1:length(Qvec))    = Qvec;
        
        % determine optimal stopping position 
        if findPick == 1 & draw < nchoices & (Qvec(1) > Qvec(3) | Qvec(2) > Qvec(3))
            pickTrial(trl)  = draw; % number of draws
            break

        elseif findPick == 1 & draw == nchoices
            pickTrial(trl)  = draw;
        end

        % 
        if draw == lseq
            d = [d; 0];
        end   
  
    end % end of draw loop
    
    % accumulate chosen urns
    [biggest_value choice(trl)] = max(Qvec);


end % end of condition trials



return