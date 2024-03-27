function io_output = run_POMDP_BeadsIo(R,thiscond_seq,thiscond_data)

% comments will go here

%% upack struct

cond        = R.cond;
totaltrials = size(thiscond_seq,2);

% what is the probability of this cond? 
if cond == 1
    thisq = R.q(1);
else 
    thisq = R.q(2);
end

R.thisq = thisq;

% what is the urntype?
urntype = thiscond_data(:,4);

% extract sequences from cell and store in matrix (26x10)
for i = 1:size(thiscond_seq,2)
    thiscond_seqmat(i,:)    = thiscond_seq{1,i};
end

% switch indexing for green urns.. (for IO)
for u = 1:length(urntype)
    if urntype(u) == 0 % if green urn switch index coding
        seq_ones = find(thiscond_seqmat(u,:) == 1);
        seq_twos = find(thiscond_seqmat(u,:) == 2);
        thiscond_seqmat(u,seq_ones) = 2;
        thiscond_seqmat(u,seq_twos) = 1;
    end 
end

% recode 2s to 0s for backward induction 
thiscond_seqmat(find(thiscond_seqmat==2)) = 0;

% run backward induction (bruno's code)
[r, Qsat]                   = backWardInduction(thiscond_seqmat, R);

% store ideal observer output
io_output.r           = r;
io_output.Qsat        = Qsat;

% loop over condition trials to compute choices, picktrials and acc
for i = 1: totaltrials
    
    choiceTrial             = find(squeeze(Qsat(i,:,3)) - max(squeeze(Qsat(i,:,1:2))') < 0); % which options this trial vec have an urn > sample
    pickTrial(i)            = choiceTrial(1); % pick the first of the choices
    [ma ma_i]               = max(squeeze(Qsat(i,pickTrial(i),:))); % which of the two urn was chosen? [based on AQ values]
    
    if (ma_i == 1 & urntype(i) == 1) | (ma_i == 2 & urntype(i) == 0)
        choice(i)           = 1;
    else
        choice(i)           = 0;
    end
end  

% get mean accuracy and draws and store in output structure
mean_choice                 = mean(choice);
io_output.accuracy          = mean_choice;
io_output.draws             = mean(pickTrial);
io_output.points            = (sum(choice==1)*R.correct) + (sum(choice==0)*R.error) + (sum(pickTrial)*R.sample);

end % end of main function 

%% run backward induction

function [reward, Qsat] = backWardInduction(thiscond_seqmat, R)

% how many trials per conditon?
Ntrials         = size(thiscond_seqmat,1);
maxDraws        = 10;

K               = 3; % choice options

reward          = zeros(Ntrials, 1);

Qsat            = zeros(Ntrials, maxDraws, K);

parfor trial = 1 : Ntrials
    
    Qsad            = zeros(maxDraws, K); % action values for this sequence will be stored here
    drawSequence    = thiscond_seqmat(trial,:);
   
%     
    for draw = 1 : maxDraws
                   
        Qsad(draw, :) = backWardUtility(drawSequence, draw, maxDraws, R);
        
    end
    
    % store this_trial Q values in a cell
    Qsat(trial,:,:) = Qsad;
    
    %%% randomize choice for symmetric values
    Qsac            = Qsad + 0.000001*randn(maxDraws, K);
    Qsa1            = Qsac(:, 1) - Qsac(:, 3);
    Qsa2            = Qsac(:, 2) - Qsac(:, 3);
    
    choice1         = find(Qsa1 > 0);
    choice2         = find(Qsa2 > 0);
    
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

end % end of backward induction function

%% run backward Utility 

function Qsa = backWardUtility(drawSequence, draw, maxDraws, R)

utility = zeros(maxDraws, maxDraws+1);

ng = sum(drawSequence(1:draw));

for drawi = maxDraws : -1 : (draw + 1)
        
    [utility] = stateUtilityBeads(utility, drawi, draw, maxDraws, ng, R);
    
end
    
Qsa = actionValueBeads(utility, R, draw, ng, draw, maxDraws);
    
end

%% run state utility beads

function utility_t = stateUtilityBeads(utility, drawi, draw, maxDraws, ng, R)

utility_t = zeros(maxDraws, maxDraws+1);

futureDraws = drawi - draw;

ndf = drawi;

for greenDraws = 0 : futureDraws
    
    ngf = ng + greenDraws;

    Qsa = actionValueBeads(utility, R, ndf, ngf, drawi, maxDraws);

    utility_t(ndf, ngf+1) = max(Qsa);        
    
end

end 

%% run action values function

function Qsa = actionValueBeads(utility, R, nd, ng, drawi, maxDraws)

pg = PG(R.thisq, nd, ng);

pb = 1 - pg;

QG = R.correct*pg + R.error*pb;
QB = R.correct*pb + R.error*pg;

if drawi < maxDraws

    QD = R.sample + pb*((1-R.thisq)*utility(nd+1, ng+1+1) +   (R.thisq)*(utility(nd+1, ng+1))) + ...
                    pg*(  (R.thisq)*utility(nd+1, ng+1+1) + (1-R.thisq)*(utility(nd+1, ng+1)));
                
else
    
    QD = 0;
    
end

Qsa = [QG; QB; QD];

end 

%% run probability function

function p = PG(q, nd, ng)

p = 1/(1 + (q/(1-q))^(nd-2*ng));

end
