function [choiceStop, choiceCont, difVal, currentRnk] = computeSecretary(Generate_params, sampleSeries, prior, N, list, Cs,  minValue)

sdevs = 8;
dx = 2*sdevs*sqrt(prior.sig)/100;
x = ((prior.mu - sdevs*sqrt(prior.sig)) + dx : dx : ...
    (prior.mu + sdevs*sqrt(prior.sig)))';

Nchoices = length(list.vals);

Nconsider = length(sampleSeries);
if Nconsider > N
    Nconsider = N;
end

difVal      = zeros(1, Nconsider);
choiceCont  = zeros(1, Nconsider);
choiceStop  = zeros(1, Nconsider);
currentRnk  = zeros(1, Nconsider);

for ts = 1 : Nconsider
    
    [expectedStop, expectedCont] = rnkBackWardInduction(sampleSeries, ts, prior, N, x, Cs,  minValue,Generate_params,list);
    % [expectedStop, expectedCont] = backWardInduction(sampleSeries, ts, priorProb, x, Cs);
    
    [rnkv, rnki] = sort(sampleSeries(1:ts), 'descend');
    z = find(rnki == ts);
    
    %     fprintf('sample %d rnk %d %.2f %.4f %.2f\n', ts, z, sampleSeries(ts), expectedStop(ts), expectedCont(ts));
    
    difVal(ts) = expectedCont(ts) - expectedStop(ts);
    
    choiceCont(ts) = expectedCont(ts);
    choiceStop(ts) = expectedStop(ts);
    
    currentRnk(ts) = z;
 
end % end of backward induction loop





retrun