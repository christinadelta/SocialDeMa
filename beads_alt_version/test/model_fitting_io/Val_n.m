function [v, d, Qvec] = Val_n(q, nd, ng, alpha, lseq, Cw, Cc, Cs)

%%% computes probability that we are drawing from green urn
pg = PG_n(q, nd, ng);
%%% probability that we are drawing from predominantly blue urn
pb = 1 - pg;

%%% cost of choosing Green, at this point
QG = Cw*pb + Cc*pg;
%%% cost of choosing blue at this point
QB = Cw*pg + Cc*pb;


if nd + 1 <= lseq
    try
        
        %%% compute value of next state given that we draw a green
        val11 = vVal_n(q, nd+1, ng+1, alpha, lseq, Cw, Cc,Cs);
        %%% compute value of next state given that we draw a blue
        val10 = vVal_n(q, nd+1, ng, alpha, lseq, Cw, Cc,Cs);
        
        %%% redundant with above
%         val00 = vVal(q, nd+1, ng+1, alpha, lseq, Cw, Cs);
%         val01 = vVal(q, nd+1, ng, alpha, lseq, Cw, Cs);

        %%% Value of action is cost to sample plus expected value of next
        %%% state
        QS = Cs + pg*(val11*q     + val10*(1-q)) +...
                  pb*(val11*(1-q) + val10*q);
              
        %%% compile action values into vector
        Qvec = [QG; QB; QS]; 
        
        if nd == 1
           fprintf(''); 
        end
    catch
        fprintf('');
    end
else
    Qvec = [QG; QB];
end

try sum(exp(alpha*Qvec));
catch; fprintf('');
end

if sum(exp(alpha*Qvec)) == 0
    fprintf('');
end

%%% softmax function to convert values to action probabilities
d = exp(alpha*Qvec)./sum(exp(alpha*Qvec));

%%% average value of this state if we take actions with probability d
v = d'*Qvec;

return