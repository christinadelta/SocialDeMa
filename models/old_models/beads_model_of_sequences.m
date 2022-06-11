function [] = beads_model_of_sequences;

sequence1 = 1101001101;
sequence2 = 1110011010;

R.correct = 0;
R.error = -1000;
R.sample = -10;
R.q = 0.6;

seq_mat = [sequence1; sequence2];

[r, Qsa] = backWardInduction(size(seq_mat,1), size(seq_mat,2), seq_mat, R);

for dri=1:size(seq_mat,1);
    choiceTrial = find(squeeze(Qsa(dri, :, 3)) - max(squeeze(Qsa(dri, :, 1:2))') < 0);   %which options this seq have an urn > sample
    pickTrial(dri) = choiceTrial(1);    %which option was the first one where saample was inferior to urn (choice)?
    [ma ma_i] = max(squeeze(Qsa(dri,pickTrial(dri),:)));    %index of max value on choice option (which urn was chosen/best)?
    choice(dri) = ma_i;  %assign chosen urn
    choice(find(choice==2)) = 0; %recode it so it can be summed
end; %seque


disp('audi5000');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% for running synthetic daata
function [reward, Qsat] = backWardInduction(Ntrials, maxDraws, drawSequence, R)

K = 3;

reward = zeros(Ntrials, 1);

Qsat = zeros(Ntrials, maxDraws, K);

for trial = 1 : Ntrials
    
    Qsad = zeros(maxDraws, 3);
    
    for draw = 1 : maxDraws
                   
        Qsad(draw, :) = backWardUtility(drawSequence(trial, :), draw, maxDraws, R);
        
    end
    
    Qsat(trial, :, :) = Qsad;
    
    %%% randomize choice for symmetric values
    Qsac = Qsad + 0.000001*randn(maxDraws, K);
    
    Qsa1 = Qsac(:, 1) - Qsac(:, 3);
    Qsa2 = Qsac(:, 2) - Qsac(:, 3);
    
    choice1 = find(Qsa1 > 0);
    choice2 = find(Qsa2 > 0);
    
    if isempty(choice1)
        choice1(maxDraws+1) = 1;
    end
    
    if isempty(choice2)
        choice2(maxDraws+1) = 1;
    end
    
    if choice1(1) < choice2(1)
        reward(trial) = 1;
    else
        reward(trial) = 0;
    end
    
            
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Qsa = backWardUtility(drawSequence, draw, maxDraws, R)

utility = zeros(maxDraws, maxDraws+1);

ng = sum(drawSequence(1:draw));

for drawi = maxDraws : -1 : (draw + 1)
        
    [utility] = stateUtilityBeads(utility, drawi, draw, maxDraws, ng, R);
    
end
    
Qsa = actionValueBeads(utility, R, draw, ng, draw, maxDraws);
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function utility_t = stateUtilityBeads(utility, drawi, draw, maxDraws, ng, R)

utility_t = zeros(maxDraws, maxDraws+1);

futureDraws = drawi - draw;

ndf = drawi;

for greenDraws = 0 : futureDraws
    
    ngf = ng + greenDraws;

    Qsa = actionValueBeads(utility, R, ndf, ngf, drawi, maxDraws);

    utility_t(ndf, ngf+1) = max(Qsa);        
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Qsa = actionValueBeads(utility, R, nd, ng, drawi, maxDraws)

pg = PG(R.q, nd, ng);

pb = 1 - pg;

QG = R.correct*pg + R.error*pb;
QB = R.correct*pb + R.error*pg;

if drawi < maxDraws

    QD = R.sample + pb*((1-R.q)*utility(nd+1, ng+1+1) +   (R.q)*(utility(nd+1, ng+1))) + ...
                    pg*(  (R.q)*utility(nd+1, ng+1+1) + (1-R.q)*(utility(nd+1, ng+1)));
                
else
    
    QD = 0;
    
end

Qsa = [QG; QB; QD];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p = PG(q, nd, ng)

p = 1/(1 + (q/(1-q))^(nd-2*ng));
