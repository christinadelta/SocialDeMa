function [] = beads_model_of_sequences;

sequence1 = [1 1 0 1 0 0 1 1 0 1];
sequence2 = [1 1 1 0 0 1 1 0 1 0];

seq_mat = [sequence1; sequence2];
% seq_mat = [sequence1];
% seq_mat = [sequence2];


 %need sequence and behavior, but behavior needs to be in its old format
           alpha = 1;  %softmax
           Cw = -1000;
           q = 0.6;
           Cs = -10;
           [ll, pickTrial, dQvec, ddec, aQvec choice] = estimateLikelihoodf(alpha,Cw,q,Cs,seq_mat,1);
           
           choice(find(choice==2)) = 0;
           all_accuracy = mean(choice==1)
           all_draws = mean(pickTrial)
           all_points = 0 + (sum(choice==0)*-1000) - sum(pickTrial)

disp('audi5000');


function [ll, pickTrial, dQvec, ddec, aQvec choice] = estimateLikelihoodf(alpha,Cw,q,Cs, sequence, aqvec_switch);

findPick = 1;

%%% length of sequence
lseq    = size(sequence, 2);

%%% number of sequences
nblocks = size(sequence, 1);

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
        if sequence(block, draw) == 1
            ng = ng + 1;
        end

        %%% alwasy increment draws
        nd = nd + 1;

        %%% compute values of each action
        [v, d, Qvec] = Val(q, nd, ng, alpha, lseq, Cw, Cs);
        
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [v] = vVal(q, nd, ng, alpha, lseq, Cw, Cs)

[v, d, Qvec] = Val(q, nd, ng, alpha, lseq, Cw, Cs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [v, d, Qvec] = Val(q, nd, ng, alpha, lseq, Cw, Cs)

%%% computes probability that we are drawing from green urn
pg = PG(q, nd, ng);
%%% probability that we are drawing from predominantly blue urn
pb = 1 - pg;

%%% cost of choosing Green, at this point
QG = Cw*pb;
%%% cost of choosing blue at this point
QB = Cw*pg;


if nd + 1 <= lseq
%     try
        
        %%% compute value of next state given that we draw a green
        val11 = vVal(q, nd+1, ng+1, alpha, lseq, Cw, Cs);
        %%% compute value of next state given that we draw a blue
        val10 = vVal(q, nd+1, ng, alpha, lseq, Cw, Cs);
        
        %%% redundant with above
%         val00 = vVal(q, nd+1, ng+1, alpha, lseq, Cw, Cs);
%         val01 = vVal(q, nd+1, ng, alpha, lseq, Cw, Cs);

        %%% Value of action is cost to sample plus expected value of next
        %%% state
        QS = Cs + pg*(val11*q     + val10*(1-q)) +...
                  pb*(val11*(1-q) + val10*q);
              
        %%% compile action values into vector
        Qvec = [QG; QB; QS]; 
%         if nd == 1
%            fprintf(''); 
%         end
%     catch
%         fprintf('');
%     end
else
    Qvec = [QG; QB];
end

% if sum(exp(alpha*Qvec)) == 0
%     fprintf('');
% end

%%% softmax function to convert values to action probabilities
d = exp(alpha*Qvec)./sum(exp(alpha*Qvec));

%%% average value of this state if we take actions with probability d
v = d'*Qvec;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p = PG(q, nd, ng)

p = 1/(1 + (q/(1-q))^(nd-2*ng));

