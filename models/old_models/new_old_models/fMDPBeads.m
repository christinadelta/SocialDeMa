function [ll, Qsad, cprob] = fMDPBeads(params, drawData, choiceData)

beta = 0.13;

if ~isempty(params)

    maxDraws = 10;
        
    Ntrials = max(size(drawData));
    
    Qsad = NaN*ones(Ntrials, maxDraws, 3);
    
    ll = 0;
    
    for trial = 1 : Ntrials
            
        trialDraws = drawData{trial}(2:end);
        
        nDraws = length(trialDraws);
        
        R.sample = params(1);
        
        if drawData{trial}(1) == 1
            R.correct = 10;
            R.error = -10;
            R.q = 0.8;
        elseif drawData{trial}(1) == 2
            R.correct = 10;
            R.error = -10;
            R.q = 0.6;
        elseif drawData{trial}(1) == 3
            R.correct = 10;
            R.error = 0;
            R.q = 0.8;
         elseif drawData{trial}(1) == 4
            R.correct = 10;
            R.error = 0;
            R.q = 0.6;
        end 

        for draw = 1 : nDraws

            Qsad(trial, draw, 1:3) = backWardUtility(trialDraws, draw, maxDraws, R)';
            
            vVec = Qsad(trial, draw, 1:3);
            
            cprob(trial, draw, :) = exp(beta*vVec)./sum(exp(beta*vVec));
            
        end
        
        seqChoice = (choiceData{trial}(end) == 2) + 1;  %%% 6 green, 2 blue;
                        
        if nDraws-1 > 0 & nDraws < maxDraws
            ll = ll - sum(log(squeeze(cprob(trial, nDraws-1, 3)))) - log(squeeze(cprob(trial, nDraws, seqChoice)));
%             ll = ll - (Qsad(trial, nDraws-1, 3) - max(Qsad(trial, nDraws-1, :))) - ...
%                       (Qsad(trial, nDraws, seqChoice) - Qsad(trial, nDraws, 3));
        elseif nDraws < maxDraws
%             ll = ll - (Qsad(trial, nDraws, seqChoice) - Qsad(trial, nDraws, 3));
            ll = ll - log(squeeze(cprob(trial, nDraws, seqChoice)));
        end  
                
    end
    
else

    Ntrials  = 4;
    maxDraws = 12;

    R.correct = 1;
    R.error   = 0;
    R.sample  = -0.025;
    R.q       = 0.6;

    drawSequence = generateDrawSequences(R.q, maxDraws, Ntrials);

    [r, Qsa] = backWardInduction(Ntrials, maxDraws, drawSequence, R);

    figure;

    for dri = 1 : 4
        subplot(2,2,dri);
        plot(1:maxDraws, squeeze(Qsa(dri, :, :)), 'LineWidth', 2);

        choiceTrial = find(squeeze(Qsa(dri, :, 3)) - max(squeeze(Qsa(dri, :, 1:2))') < 0);
        text(choiceTrial(1), 0.1, 'X', 'horizontalalignment', 'center');

    %     legend('Green', 'Blue', 'Draw');

        hold on;
        plot(1:maxDraws, drawSequence(dri, :)*0.8 + 0.1, 'ro');

        axis([0 maxDraws + 1 0 1]);

    end

    figure;

    R.sample = -0.005;

    [r, Qsa] = backWardInduction(Ntrials, maxDraws, drawSequence, R);

    for dri = 1 : 4
        subplot(2,2,dri);
        plot(1:maxDraws, squeeze(Qsa(dri, :, :)), 'LineWidth', 2);
    %     legend('Green', 'Blue', 'Draw');

        choiceTrial = find(squeeze(Qsa(dri, :, 3)) - max(squeeze(Qsa(dri, :, 1:2))') < 0);
        text(choiceTrial(1), 0.1, 'X', 'horizontalalignment', 'center');


        hold on;
        plot(1:maxDraws, drawSequence(dri, :)*0.8 + 0.1, 'ro');

        axis([0 maxDraws + 1 0 1]);

    end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function drawSequence = generateDrawSequences(q, maxDraws, Ntrials)

drawSequence = (rand(Ntrials, maxDraws) < q);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% for running synthetic daata
function [reward, Qsat] = backWardInduction(Ntrials, maxDraws, drawSequence, R)

K = 3;

reward = zeros(Ntrials, 1);

Qsat = zeros(Ntrials, maxDraws, K);

parfor trial = 1 : Ntrials
    
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