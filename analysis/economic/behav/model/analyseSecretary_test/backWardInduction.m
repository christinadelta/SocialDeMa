function [expectedStop, expectedCont, expectedUtility] = backWardInduction(sampleSeries, ts, prior, x, Cs)

N = length(sampleSeries);

data.n  = ts;

% if ts > 0
data.sig = var(sampleSeries(1:ts));
data.mu = mean(sampleSeries(1:ts));
% else
%     data.sig = priorProb.sig;
%     data.mu  = priorProb.mu;
% end

utStop  = x;

utCont = zeros(N, 1);
utility = zeros(length(x), N);

if ts == 0
    ts = 1;
end

for ti = N : -1 : ts
    
    expData = data;
    expData.n = ti;
    
    [postProb] = normInvChi(prior, expData);
    
    px = posteriorPredictive(x, postProb);
    
    px = px/sum(px);
    
    if ti == N
        utCont(ti) = -Inf;
    else
        utCont(ti) = sum(px.*utility(:, ti+1)) - Cs;
    end
    
    if ti == ts
        utility(:, ti)   = ones(size(utStop, 1), 1)*max([sampleSeries(ts) utCont(ti)]);
        expectedStop(ti) = sampleSeries(ts);
    else
        utility(:, ti)   = max([utStop ones(size(utStop, 1), 1)*utCont(ti)], [], 2);
        expectedStop(ti) = sum(px.*utStop);
    end
    
    expectedUtility(ti) = sum(px.*utility(:,ti));
    expectedCont(ti)    = utCont(ti);
    
end

return