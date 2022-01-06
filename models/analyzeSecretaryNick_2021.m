function [choiceStop, choiceCont, difVal] = analyzeSecretary_Nick2021(Generate_params,list);
% dataPrior, list, input_params, distOptions,minValue,param_to_fit,experiment, subject)

%2021:updated for use witah Param_rover_v1.m
%Written assuming it operating on exactly onse sequence at a time

%Assign params. Hopefully everything about model is now pre-specified
%(I would update mean and variance directly from Generate_params 
prior.mu    = Generate_params.PriorMean + Generate_params.model(Generate_params.current_model).BP; %prior mean offset is zero unless biased prior model
prior.sig   = Generate_params.PriorVar + Generate_params.model(Generate_params.current_model).BPV;                        %Would a biased variance model be a distinct model?
if prior.sig < 1; prior.sig = 1; end;   %It can happen randomly that a subject has a low variance and subtracting the bias gives a negative variance. Here, variance is set to minimal possible value.
prior.kappa = Generate_params.model(Generate_params.current_model).kappa;   %prior mean update parameter
prior.nu    = Generate_params.model(Generate_params.current_model).nu;

%Cost to sample
Cs = Generate_params.model(Generate_params.current_model).Cs;      %Will already be set to zero unless Cs model

%Bin edges
minValue =  Generate_params.binEdges_reward;

%If there are BV parameters specified, then warp the inputs
% if ~isnan(Generate_params.model(Generate_params.current_model).BVslope);
if Generate_params.model(Generate_params.current_model).identifier == 4;   %If identifier is BV
    list.allVals = ...
        (Generate_params.BVrange(2) - Generate_params.BVrange(1)) ./ ...
        (1+exp(-Generate_params.model(Generate_params.current_model).BVslope*(list.allVals-Generate_params.model(Generate_params.current_model).BVmid))); %do logistic transform
end;

list.vals = list.allVals';

% if list.flip == -1
%     sampleSeries = -(dataList - mean(dataList)) + mean(dataList);
% else
    sampleSeries = list.vals;
% end

N = Generate_params.seq_length;

% list.param_to_fit = param_to_fit;

[choiceStop, choiceCont, difVal, currentRnk] = computeSecretary(Generate_params, sampleSeries, prior, N, list, Cs,minValue);
% 
% if list.optimize == 1
%     z = find(difVal < 0);
%     [~, rnki] = sort(sampleSeries, 'descend');
%     rnkValue = find(rnki == z(1));
%     
%     winnings = (rnkValue == 1)*5 + (rnkValue == 2)*2 + (rnkValue == 3)*1;
% else
%     winnings = 0;
%     rnkValue = -1*ones(length(list.vals), 1);
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [choiceStop, choiceCont, difVal, currentRnk] = computeSecretary(Generate_params, sampleSeries, priorProb, N, list, Cs,  minValue)

sdevs = 8;
dx = 2*sdevs*sqrt(priorProb.sig)/100;
x = ((priorProb.mu - sdevs*sqrt(priorProb.sig)) + dx : dx : ...
    (priorProb.mu + sdevs*sqrt(priorProb.sig)))';

Nchoices = length(list.vals);

% if list.optimize == 1
%     Nconsider = length(list.allVals);
% else
    Nconsider = length(sampleSeries);
    if Nconsider > N
        Nconsider = N;
    end
% end

difVal = zeros(1, Nconsider);
choiceCont = zeros(1, Nconsider);
choiceStop = zeros(1, Nconsider);
currentRnk = zeros(1, Nconsider);

for ts = 1 : Nconsider
    
    [expectedStop, expectedCont] = rnkBackWardInduction(sampleSeries, ts, priorProb, N, x, Cs,  minValue,Generate_params,list);
    %     [expectedStop, expectedCont] = backWardInduction(sampleSeries, ts, priorProb, x, Cs);
    
    [rnkv, rnki] = sort(sampleSeries(1:ts), 'descend');
    z = find(rnki == ts);
    
    %     fprintf('sample %d rnk %d %.2f %.4f %.2f\n', ts, z, sampleSeries(ts), expectedStop(ts), expectedCont(ts));
    
    difVal(ts) = expectedCont(ts) - expectedStop(ts);
    
    choiceCont(ts) = expectedCont(ts);
    choiceStop(ts) = expectedStop(ts);
    
    currentRnk(ts) = z;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [expectedStop, expectedCont, expectedUtility] = backWardInduction(sampleSeries, ts, priorProb, x, Cs)

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
    
    [postProb] = normInvChi(priorProb, expData);
    
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [expectedStop, expectedCont, expectedUtility] = rnkBackWardInduction(sampleSeries, ts, priorProb, ...
    listLength, x, Cs,  minValue,Generate_params,list)

% maxPayRank = 3;
% payoff = [5 3 1 0 0 0 0 0 0 0 0 0 0 0 0 0 ];
% % payoff = [1 0 0 0 0 0];



N = listLength;
Nx = length(x);

% payoff = sort(sampleSeries,'descend')';
% payoff = [N:-1:1];
% payoff = (payoff-1)/(N-1);

% if list.param_to_fit == 9;  %if BR model

if Generate_params.model(Generate_params.current_model).identifier == 5;    
    %     Generate_params.model(Generate_params.current_model).BR.slope = NaN;
%     Generate_params.model(Generate_params.current_model).BR.mid = NaN;
    rewards = ...
        (Generate_params.BVrange(2) - Generate_params.BVrange(1)) ./ ...
        (1+exp(-Generate_params.model(Generate_params.current_model).BRslope*(list.allVals-Generate_params.model(Generate_params.current_model).BRmid))); %do logistic transform

%     rewards = (params(1)-params(2))./(1+exp(-params(3)*(sampleSeries-params(4)))); %do logistic transform
    payoff = sort(rewards,'descend');
    payoff = (payoff-1)/(Generate_params.BVrange(2)-1);  %normalises bin values to 0 to 1 scale (for some reason, presumable so I can compare models with different reward functions by placing them on same scale, although I'm not ure scale matters. Will keep for now)
    
else;
    
    % % %bins
    temp = sort(sampleSeries,'descend')';   %sort the sample values
    [dummy,payoff] = histc(temp, [minValue(1:end-1) Inf]);  %Slot the sequence values into the 6 (or howevermany) bins
    nbins = size(minValue,2)-1;     %recover number of bins if I didn't pass it through
    payoff = (payoff-1)/(nbins-1);  %normalises bin values to 0 to 1 scale (for some reason, presumable so I can compare models with different reward functions by placing them on same scale, although I'm not ure scale matters. Will keep for now)
    
end;

% %normalised rating value
% payoff = sort(sampleSeries,'descend')'; %assign actual values to payoff
% payoff = (payoff-0)/(minValue(end) - 0);    %normalise seq values between zero and 1 relative to maximally rated face




% payoff(find(payoff~=8))=0.0000000000000000000000000000000001;
% payoff(find(payoff==8))=100000000000000000000000000000;
%
% %bound values between zero and 1
% if numel(minValue) > 2;
% payoff = ((payoff-0)/((numel(minValue)-1)-0));
% end;
% payoff = payoff.^40;

maxPayRank = numel(payoff);
payoff = [payoff zeros(1, 20)];



data.n  = ts;

% if ts > 0

% if distOptions == 0
    data.sig = var(sampleSeries(1:ts));
    data.mu = mean(sampleSeries(1:ts));
    
% else
%     data.mu = mean(sampleSeries(1:ts));
%     data.sig = data.mu;
% end
% else
%     data.sig = priorProb.sig;
%     data.mu  = priorProb.mu;
% end

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

[postProb] = normInvChi(priorProb, data);

% if list.param_to_fit == 11; %if optimism model .. 
%     postProb.mu = postProb.mu+params(1); %...Then add constant to the posterior mean
    postProb.mu ...
        = postProb.mu ...
        + Generate_params.model(Generate_params.current_model).optimism; %...Then add constant to the posterior mean (will be zero if not optimism model)
% end;

px = posteriorPredictive(x, postProb);
px = px/sum(px);

Fpx = cumsum(px);
cFpx = 1 - Fpx;

for ti = N : -1 : ts
    
    if ti == N
        utCont = -Inf*ones(Nx, 1);
    elseif ti == ts
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
    
    %     psi(:,1) = (Fpx.^(td));
    %     psi(:,2) = td*(Fpx.^(td-1)).*(cFpx);
    %     psi(:,3) = (td*(td-1)/2)*(Fpx.^(td-2)).*(cFpx.^2);
    
    for ri = 1 : maxPayRank+1
        
        z = find(x < rnkv(ri) & x >= rnkv(ri+1));
        utStop(z) = ps(z, 1:maxPayRank)*(payoff(1+(ri-1):maxPayRank+(ri-1))');
        
    end
    
    if sum(isnan(utStop)) > 0
        fprintf('Nan in utStop');
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
    
    utCont = utCont - Cs;
    
    utility(:, ti)      = max([utStop utCont], [], 2);
    expectedUtility(ti) = px'*utility(:,ti);
    
    expectedStop(ti)    = px'*utStop;
    expectedCont(ti)    = px'*utCont;
    
    %     subplot(2,1,1);
    %     plot(x, utStop, x, utCont, x, utility(:, ti));
    %
    %     subplot(2,1,2);
    %     plot(x, Fpx);
    %
    %     fprintf('');
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function utCont = computeContinue(utility, postProb0, x, ti)

postProb0.nu = ti-1;

utCont = zeros(length(x), 1);

% pspx = zeros(length(x), length(x));

expData.n   = 1;
expData.sig = 0;

for xi = 1 : length(x)
    
    expData.mu  = x(xi);
    
    postProb = normInvChi(postProb0, expData);
    spx = posteriorPredictive(x, postProb);
    spx = (spx/sum(spx));
    
    %     pspx(:, xi) = spx;
    
    utCont(xi) = spx'*utility;
    
end

% subplot(2,2,1);
% plot(x, pspx(:, 1:100:end));
%
% subplot(2,2,2);
% plot(x, utility);
% title(ti);
%
% subplot(2,2,3);
% plot(x, utCont);
% axis([min(x) max(x) 0.5 3]);
%
% fprintf('');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [postProb] = normInvChi(prior, data)

postProb.nu    = prior.nu + data.n;

postProb.kappa = prior.kappa + data.n;

postProb.mu    = (prior.kappa/postProb.kappa)*prior.mu + (data.n/postProb.kappa)*data.mu;

postProb.sig   = (prior.nu*prior.sig + (data.n-1)*data.sig + ...
    ((prior.kappa*data.n)/(postProb.kappa))*(data.mu - prior.mu).^2)/postProb.nu;

if data.n == 0
    postProb.sig = prior.sig;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function prob_y = posteriorPredictive(y, postProb)

tvar = (1 + postProb.kappa)*postProb.sig/postProb.kappa;

sy = (y - postProb.mu)./sqrt(tvar);

prob_y = tpdf(sy, postProb.nu);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [z, fz] = integrateSampling(x, postProb)

sdevs = 8;

maxx = postProb.mu + sdevs*sqrt(postProb.sig);
minx = postProb.mu - sdevs*sqrt(postProb.sig);

dx = (maxx - minx)/10000;

xv = ((minx + dx) : dx : maxx)';

pv = posteriorPredictive(xv, postProb);

pv = pv/sum(pv);

[v, vi] = min(abs(x - xv));

z = sum(pv((vi+1):end));

fz = sum(xv((vi+1):end).*pv(((vi+1):end)))/z;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function probLarger = forwardInduction(sampleSeries, priorProb)

N = length(sampleSeries);

for t = 1 : N
    
    utStop = sampleSeries(t);
    
    data.mu = mean(sampleSeries(1:t));
    data.n  = t;
    
    if t > 1
        data.sig = var(sampleSeries(1:t));
    else
        data.sig = 0;
    end
    
    [postProb] = normInvChi(priorProb, data);
    
    probLarger(t) = 1 - (1 - integrateSampling(sampleSeries(t), postProb)).^(N-t);
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function extraCode()

% %%% Need to see if this is a candidate...
% [~, rnki] = sort(sampleSeries(1:ts), 'descend');
% rnk = rnki(end);
%
%
%
% if rnk > 3
%     expectedReward = 0;
% else
%
%     [pgreater, ~] = integrateSampling(sampleSeries(ts), postProb);
%     pless = 1 - pgreater;
%
%     rd = N - ts;
%
%     ps = zeros(length(payoff), 1);
%
%     ps(1+(rnk-1)) = pless^(rd);
%     ps(2+(rnk-1)) = rd*(pless^(rd-1))*(pgreater);
%     ps(3+(rnk-1)) = (rd*(rd-1)/2)*(pless^(rd-2))*(pgreater^2);
%
%     expectedReward = payoff*ps;
%
% end








