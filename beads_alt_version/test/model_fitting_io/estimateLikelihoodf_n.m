function [ll, pickTrial, dQvec, ddec, aQvec] = estimateLikelihoodf_n(params, seq_mat, choiceVec, fixedParams, findPick)

%%% extract free parameters 
% Cw = params(1);
Cs      = params(1);

% extract fixed parameters
alpha   = fixedParams(1);
q       = fixedParams(2);
Cw      = fixedParams(3);
Cc      = fixedParams(4);
cond    = fixedParams(5);


%%% number of sequences
nblocks = size(choiceVec, 2);

%%% intialize log likelihood to zero
ll = 0;

for block = 1 : nblocks % blocks are the number of sequences per condition (26)
    
    % extract this sequence
    this_seq = seq_mat(block,:);
    
    %%% choices of subject for this sequence of draws
    thischoice = choiceVec{cond,block};
    
    %%% length of sequence
    lseq    = size(this_seq, 2);

    %%% number of choices for this sequence
    nchoices = size(thischoice, 1);
    
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
        [v, d, Qvec] = Val_n(q, nd, ng, alpha, lseq, Cw, Cc, Cs);
        
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
        
        if thischoice(draw,:)*d == 0;
            fprintf('missing data sequence %d   ', block);
            thischoice(draw,:) = [1/3 1/3 1/3];
        end
            

        %%% update log likelihood
        try
            ll = ll - log(thischoice(draw,:)*d);
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