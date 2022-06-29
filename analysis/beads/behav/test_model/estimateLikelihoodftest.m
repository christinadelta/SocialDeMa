function [ll, pickTrial, dQvec, ddec, aQvec choice] = estimateLikelihoodftest(alpha,Cw,q,Cs, sequence, aqvec_switch);

findPick = 1;

%%% length of sequence
lseq    = 10;

%%% number of sequences
nblocks = size(sequence, 2); 

%%% intialize log likelihood to zero
ll = 0;

for block = 1 : nblocks
    
%     %% choices of subject for this sequence of draws
%     choiceVec = setData{block};

    %%% number of choices for this sequence
%     nchoices = size(choiceVec, 1);
    nchoices = lseq;
    
    %%% initially nd (draws) == 0 and ng (green marbles) == 0
    ng = 0;
    nd = 0;
    
    %%% Qvec is values of each action
    dQvec = [];
    
    %%% corresponding probabilities generated with softmax and alpha
    ddec  = [];

    %%% loop over draws for this sequence of draws
    for draw = 1 : nchoices

        %%% check if we got a green or blue marble
        if sequence{block}(draw) == 1
            ng = ng + 1;
        end

        %%% alwasy increment draws
        nd = nd + 1;

        %%% compute values of each action
        [v, d, Qvec] = Valtest(q, nd, ng, alpha, lseq, Cw, Cs);
        
        %%% keep track of values across sequnce of draws
        dQvec(draw, 1:length(Qvec)) = Qvec;
        
        %%% keep track of choice probabilities across sequence 
        ddec (draw, 1:length(d))    = d;
        
        %%% Nick, add this line to your code and also return it.
        aQvec{block}(draw, 1:length(Qvec)) = Qvec;
        
        %%% if trying to determine optimal stopping position (findpick ==
        %%% 1) then see if we should stop
        if findPick == 1 & draw < nchoices & (Qvec(1) > Qvec(3) | Qvec(2) > Qvec(3))
            pickTrial(block) = draw;
            break
        elseif findPick == 1 & draw == nchoices
            pickTrial(block) = draw;
        end

        %%% if at end of sequence d is two element and have to tack on for
        %%% syntax purposes
        if draw == lseq
            d = [d; 0];
        end

        %%% update log likelihood
%         try
%             ll = ll - log(choiceVec(draw,:)*d);
%         catch
%             fprintf('');
%         end
        
%         if findPick == 0
% 
%             if length(Qvec) == 3
%                 fprintf('%d\t', block);
%                 for di = 1 : 3
%                     fprintf('%.2f/%.2f\t', Qvec(di), d(di));
%                 end
%                 fprintf('%d pg %.2f\n', sequence(block, draw), PG(q, nd, ng));
%             end
%             
%         end

    end;    %loop through draws
    
    %accumlate chosen urns
    [biggest_value choice(block)] = max(Qvec);
    
    
%     subplot(2,2,block);
%     plot(dQvec);
%     legend('G', 'B', 'S');
    
end;    %loop through sequences

% if findPick == 0
%     fprintf('ll %.2f Cw %.1f Cs %.1f alpha %.2f\n', ll, Cw, Cs, alpha);
% end



return