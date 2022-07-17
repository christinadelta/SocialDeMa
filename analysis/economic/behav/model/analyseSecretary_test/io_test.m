% TEST IDEAL OBSERVER - BEST CHOICE TASK VERSION

% this version is a test to see if the ideal observer works 
% this operates on one sequence at a time

% the parameters that the ideal observer should need are:
% 1. sequences
% 2. mean and variance of raw prices 
% 3. mean and variance of ratings 
% 4. kappa=2 and nu=1 (dof)
% 5. Cs (set to 0)?
% 6. ranks?
% 7. vector with (within sequence) responses?

% Assign params. Hopefully everything about model is now pre-specified
% (I would update mean and variance directly from Generate_params 
prior.mu    = Generate_params.PriorMean + Generate_params.model(Generate_params.current_model).BP; % prior mean offset is zero unless biased prior model
prior.sig   = Generate_params.PriorVar + Generate_params.model(Generate_params.current_model).BPV; % Would a biased variance model be a distinct model?
if prior.sig < 1; prior.sig = 1; end % It can happen randomly that a subject has a low variance and subtracting the bias gives a negative variance. Here, variance is set to minimal possible value.
prior.kappa = Generate_params.model(Generate_params.current_model).kappa; % prior mean update parameter
prior.nu    = Generate_params.model(Generate_params.current_model).nu;

% Cost to sample
Cs = Generate_params.model(Generate_params.current_model).Cs; % Will already be set to zero unless Cs model

% Bin edges
minValue =  Generate_params.binEdges_reward;

% If there are BV parameters specified, then warp the inputs
% if ~isnan(Generate_params.model(Generate_params.current_model).BVslope);
if Generate_params.model(Generate_params.current_model).identifier == 4;   %If identifier is BV
    list.allVals = ...
        (Generate_params.BVrange(2) - Generate_params.BVrange(1)) ./ ...
        (1+exp(-Generate_params.model(Generate_params.current_model).BVslope*(list.allVals-Generate_params.model(Generate_params.current_model).BVmid))); %do logistic transform
end

list.vals = list.allVals';
sampleSeries = list.vals;
N = Generate_params.seq_length;

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











