function [ll, pickTrial, dQvec, ddec, aQvec] = estimateLikelihoodf(params, sequence, setData, fixedParams, findPick, urntype)

%%% extract free parameters 
% Cw = params(1);
Cs      = params(1);
% comment this out when running through beads preprocessing script
% setData = choiceVec;

% extract fixed parameters
alpha   = fixedParams(1);
q       = fixedParams(2);
Cw      = fixedParams(3);
Cc      = fixedParams(4);
Cd      = fixedParams(5);
cond    = fixedParams(6);
% Cs = fixedParams(3);

%%% number of sequences
nblocks = size(setData, 2);

%%% intialize log likelihood to zero
ll = 0;

for block = 1 : nblocks % blocks are the number of sequences per condition (26)
    
    % extract this sequence
    this_seq = sequence{block};
    
    %%% choices of subject for this sequence of draws
    choiceVec = setData{block};
    
    %%% length of sequence
    lseq    = size(this_seq, 2);

    %%% number of choices for this sequence
    nchoices = size(choiceVec, 1);
    
    %%% NEED TO ASK NICK ABOUT THIS
    % now, i'm not sure if this part is needed but, if this is a green type
    % sequence, swipe the codes (in the sequence)
    if urntype(block) == 0
        seq_ones            = find(this_seq == 1);
        seq_twos            = find(this_seq == 2);
        this_seq(seq_ones)  = 2;
        this_seq(seq_twos)  = 1;
    end

    %%% initially nd (draws) == 0 and ng (green marbles) == 0
    ng = 0;
    nd = 0;
    
    %%% Qvec is values of each action
    dQvec = [];
    
    %%% corresponding probabilities generated with softmax and alpha
    ddec  = []; 
        
%     if isempty(find( mean( choiceVec') == 0 )) == 0; 
%         disp(sprintf('missing response skipping sequence %d for this type', block)); continue; 
%     end;

    %%% loop over draws for this sequence of draws
    for draw = 1 : nchoices

        %%% check if we got a green or blue marble
        if this_seq(draw) == 1
            ng = ng + 1;
        end

        %%% alwasy increment draws
        nd = nd + 1;

        %%% compute values of each action
        [v, d, Qvec] = Val(q, nd, ng, alpha, lseq, Cw, Cc,Cd, Cs);
        
        %%% keep track of values across sequnce of draws
        dQvec(draw, 1:length(Qvec)) = Qvec;
        
        %%% keep track of choice probabilities across sequence 
        ddec (draw, 1:length(d))    = d;
        
        %%% Nick, add this line to your code and also return it.
        aQvec{block}(draw, 1:length(Qvec)) = Qvec;
        
        %%% if trying to determine optimal stopping position (findpick ==
        %%% 1) then see if we should stop
        if findPick == 1 && draw < nchoices && (Qvec(1) > Qvec(3) || Qvec(2) > Qvec(3))
            pickTrial(block) = draw;
            break
        elseif findPick == 1 && draw == nchoices
            pickTrial(block) = draw;
        end

        %%% if at end of sequence d is two element and have to tack on for
        %%% syntax purposes
        if draw == lseq
            d = [d; 0];
        end
        
        if choiceVec(draw,:)*d == 0;
            fprintf('missing data sequence %d   ', block);
            choiceVec(draw,:) = [1/3 1/3 1/3];
        end
            

        %%% update log likelihood
        try
            ll = ll - log(choiceVec(draw,:)*d);
        catch
            fprintf('');
        end
        
        % fprintf( 'draw: %d ll: %.2f', nd, ll);
        
        % save all ll's to take a look at them 
%         all_ll{1,block}(draw,1) = ll;
        

    end
    
%     subplot(2,2,block);
%     plot(dQvec);
%     legend('G', 'B', 'S');

%     subplot(6,4,block);
%     h = plot(1:size(dQvec,1), dQvec);
%     legend('B', 'G', 'D');
%     set(gca, 'Fontname', 'Ariel', 'FontSize', 6);
%     set(h, 'MarkerSize',6, 'marker', 'o');
%     
end

if findPick == 0
    pickTrial = [];
    fprintf('ll %.2f Cw %.1f p %.2f Cs %.2f alpha %.2f\n', ll, Cw, Cc, q, Cs, alpha);
end

return