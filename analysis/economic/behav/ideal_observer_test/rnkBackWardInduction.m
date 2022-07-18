function [expectedStop, expectedCont, expectedUtility] = rnkBackWardInduction(sampleSeries, ts, prior, ...
    N, x, Cs,  minValue,Generate_params,list)

% maxPayRank = 3;
% payoff = [5 3 1 0 0 0 0 0 0 0 0 0 0 0 0 0 ];
% % payoff = [1 0 0 0 0 0];



% N = listLength;
Nx = length(x);

payoff = sort(sampleSeries,'descend')';
payoff = [N:-1:1];
payoff = (payoff-1)/(N-1);

% % %bins
% temp = sort(sampleSeries,'descend')';   % sort the sample values
% [dummy,payoff] = histc(temp, [minValue(1:end-1) Inf]); % Slot the sequence values into the 6 (or howevermany) bins
% nbins = size(minValue,2)-1; % recover number of bins if I didn't pass it through
% payoff = (payoff-1)/(nbins-1); % normalises bin values to 0 to 1 scale (for some reason, presumable so I can compare models with different reward functions by placing them on same scale, although I'm not ure scale matters. Will keep for now)
% 

% %normalised rating value
% payoff = sort(sampleSeries,'descend')'; %assign actual values to payoff
% payoff = (payoff-0)/(minValue(end) - 0);    %normalise seq values between zero and 1 relative to maximally rated face

% payoff(find(payoff~=8))=0.0000000000000000000000000000000001;
% payoff(find(payoff==8))=100000000000000000000000000000;
%
% % bound values between zero and 1
% if numel(minValue) > 2;
% payoff = ((payoff-0)/((numel(minValue)-1)-0));
% end;
% payoff = payoff.^40;

maxPayRank  = numel(payoff);
payoff      = [payoff zeros(1, 20)];
data.n      = ts;

% this part takes the mean and variance of raw prices in a sequence and not
% of all raw prices 
data.sig    = var(sampleSeries(1:ts));
data.mu     = mean(sampleSeries(1:ts));

utCont  = zeros(length(x), 1);
utility = zeros(length(x), N);

if ts == 0
    ts = 1;
end

[rnkvl, rnki] = sort(sampleSeries(1:ts), 'descend');
z = find(rnki == ts);
rnki = z;

ties = 0;
if length(unique(sampleSeries(1:ts))) < ts
    ties = 1;
end

mxv = ts;
if mxv > maxPayRank
    mxv = maxPayRank;
end

% rnkv = [Inf*ones(1,1); rnkvl(1:mxv)'; -Inf*ones(20, 1)];
rnkv = [Inf*ones(1,1); rnkvl(1:mxv); -Inf*ones(20, 1)];

[postProb] = normInvChi(prior, data);

% for ideal observer - this is not needed
% postProb.mu ...
%         = postProb.mu ...
%         + Generate_params.model(Generate_params.current_model).optimism; %...Then add constant to the posterior mean (will be zero if not optimism model)

% I am not sure if this is needed for the ideal observer
postProb.mu ...
    = postProb.mu ...
    + 0; %...Then add constant to the posterior mean (will be zero if not optimism model)
    


px = posteriorPredictive(x, postProb);
px = px/sum(px);

% compute cumulative sum
Fpx = cumsum(px);
cFpx = 1 - Fpx;

for ti = N : -1 : ts
    
    if ti == N % if this is the last draw
        utCont = -Inf*ones(Nx, 1);
    elseif ti == ts % if this is the first draw
        utCont = ones(Nx, 1)*sum(px.*utility(:, ti+1));
    else
        utCont = computeContinue(utility(:, ti+1), postProb, x, ti);
    end
    
    %%%% utility when rewarded for best 3, $5, $2, $1
    utStop = NaN*ones(Nx, 1);
    
    rd = N - ti; %%% remaining draws
    id = max([(ti - ts - 1) 0]); %%% intervening draws
    td = rd + id;
    ps = zeros(Nx, maxPayRank);
    
    for rk = 0 : maxPayRank-1
        
        pf = prod(td:-1:(td-(rk-1)))/factorial(rk);
        
        ps(:, rk+1) = pf*(Fpx.^(td-rk)).*((cFpx).^rk);
        
    end
    
    for ri = 1 : maxPayRank+1
        
        z = find(x < rnkv(ri) & x >= rnkv(ri+1));
        utStop(z) = ps(z, 1:maxPayRank)*(payoff(1+(ri-1):maxPayRank+(ri-1))');
        
    end
    
    if ti == ts
        [zv, zi] = min(abs(x - sampleSeries(ts)));
        if zi + 1 > length(utStop)
            %             fprintf('accessing utStop at %d value x %.2f\n', zi, x);
            zi = length(utStop) - 1;
        end
        
        %         if rnki > 3 & utStop(zi+1) > 0.0001 & ties == 0
        % %             fprintf('expectedReward %.9f\n', utStop(zi+1));
        %         end
        
        utStop = utStop(zi+1)*ones(Nx, 1);
        
    end
    
    utCont              = utCont - Cs;
    
    utility(:, ti)      = max([utStop utCont], [], 2);
    expectedUtility(ti) = px'*utility(:,ti);
    
    expectedStop(ti)    = px'*utStop;
    % expectedCont(:,ti)    = px'*utCont;
    expectedCont(ti)    = px'*utCont;
    
    
end % end of ti loop

    
    
return 