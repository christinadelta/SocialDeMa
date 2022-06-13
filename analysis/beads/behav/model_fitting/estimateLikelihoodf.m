function [ll, picktrial, dQvec, ddec, aQvec] = estimateLikelihoodf(params, sequence, setdata, fixedparams, findpick)

%%% extract parameters from parameter vector
Cw = params(1);
% Cs = params(2);
% alpha  = params(3);
% q = fixedparams;

alpha   = fixedparams(1);
q       = fixedparams(2);
Cs      = fixedparams(3);

lseq    = size(sequence, 2);    % length of sequence
ll      = 0;                    % init log likelihood to 0

% rename current sequence
this_sequence       = sequence;

% choices of subject for current sequence of draws 
choiceVec           = setdata;

% number of choices for this sequence (number of draws)
nchoices            = size(choiceVec,1);

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
    [v, d, Qvec]    = Val(q, numDraws, numGreen, alpha, lseq, Cw, Cs);

    aQvec(j, 1:length(Qvec)) = Qvec;
    
    % determine optimal stopping position 
    if findpick == 1 & j < nchoices & (Qvec(1) > Qvec(3) | Qvec(2) > Qvec(3))
        picktrial   = j; % number of draws
        break

    elseif findPick == 1 & j == nchoices
        picktrial = j;
    end

    % 
    if j == lseq
        d = [d; 0];
    end  
    
    % I need to work on this
%     if choiceVec(j,:) * d == 0
%         fprintf('missing data sequence')
%         coiceVec(j,:) = [1/3 1/3 1/3];
%     end

    % update log likelihood
    try 
        ll = ll - log(choiceVec(j,:) * d);
    catch
        fprintf('');
    end

    
end % end of draws for loop

% if findPick == 0
%     fprintf('ll %.2f Cw %.1f Cs %.1f alpha %.2f\n', ll, Cw, Cs, alpha);
% end


return