%% RUN IDEAL OBSERVER %%

% this script tests the ideal observer with the economic task. Here I
% define all the params needed to run the io model. Once that works, I will
% implement it through the main preprocessing script 
% The script operates on one participant at a time. 

% The input/data needed are taken from economic_prepro_v2.m (all data but
% extraction of in-sequences info -- last analysis part) 

% initiate a few parameters
log_or_not                  = 1; % log or not?
subjects                    = 1; % one subject at a time
subI                        = 1; % test subject

%% EXTRACT SUBJECT DATA AND SAVE IN PARAMS %%

% for now we need to extract this_subject phase 1 ratings [400x1], sequences
% [1x40x10] and number of samples [40x1].  
subrates                    = allsubs_ratings{1,subI};
subsequences                = allsubs_price_sequences{1,subI};
subsamples                  = allsubs_data{1,subI}.samples;
subranks                    = allsubs_data{1,subI}.rank;
subratesequences            = allsubs_rate_sequences{1,subI};
subprices                   = allsubs_prices{1,subI};

% 1. get ranks? this seems a bit confusing, because I already have the
% rank of options that each subject chose in allsubs_data.rank struct
% field. Also I don't get the code in line 93 in Nick's code
% for i = 1:length(subsequences)
%     sequenceranks{1,i} = tiedrank(subsequences{1,i})';
% end
% 
% 2 log_or_not 
if log_or_not == 1
    
    Generate_params.ratings(:,subI)             = log(subrates(:,1));
    Generate_params.prices(:,subI)              = log(subprices(:,2)); % for computing prior mean & variance

    for i = 1:length(subsequences)
        Generate_params.seq_vals(i,:,subI)      = log(subsequences{1,i})';
        Generate_params.rateseq_vals(i,:,subI)  = log(subratesequences{1,i})';
    end
else
    Generate_params.ratings(:,subI)             = subrates;
    Generate_params.prices(:,subI)              = (subprices(:,2));

    for i = 1:length(subsequences)
        Generate_params.seq_vals(i,:,subI)      = subsequences{1,i}';
        Generate_params.rateseq_vals(i,:,subI)  = (subratesequences{1,i})';
    end
end

Generate_params.num_samples(:,subI)             = subsamples(:,1);
% how do I save ranks here? Is it the rank of the option that the
% participant chose? 
Generate_params.ranks(:,sub)                    = subranks; % ranks may not be needed

% define a few more params 
Generate_params.num_subs                    = size(Generate_params.seq_vals,3);
Generate_params.num_seqs                    = size(Generate_params.seq_vals,1);
Generate_params.seq_length                  = size(Generate_params.seq_vals,2);
Generate_params.num_vals                    = size(Generate_params.ratings,1);
% Generate_params.num_vals                    = size(Generate_params.prices,1);
Generate_params.rating_bounds               = [1 100]; % What is min and max of rating scale? (Works for big trust anyway)
if log_or_not == 1
    Generate_params.rating_bounds           = log(Generate_params.rating_bounds);
end

Generate_params.BVrange                     = Generate_params.rating_bounds;
Generate_params.nbins_reward                = numel(Generate_params.rating_bounds(1):Generate_params.rating_bounds(2));  %This should effectuvely remove the binning
Generate_params.binEdges_reward             = ...
    linspace(...
    Generate_params.BVrange(1) ...
    ,Generate_params.BVrange(2)...
    ,Generate_params.nbins_reward+1 ...
    ); % organise bins by min and max

Generate_params.kappa                       = 2; 
Generate_params.nu                          = 1;
Generate_params.Cs                          = 0;

%% Generate model data %% 

% I think I will add this to a function when formally coding the io model 

% model operates one one sequence at a time
for sequence = 1:Generate_params.num_seqs
    
    % squeezing here if more than one subjects
    list.allVals                = squeeze(Generate_params.seq_vals(sequence,:,subI)); % squeezing may not be needed
    
    % should get mean and variance of sequence or whole dataset?
    Generate_params.PriorMean   = mean(Generate_params.prices(:,subI));
    Generate_params.PriorVar    = var(Generate_params.prices(:,subI));
%     Generate_params.PriorMean   = mean(Generate_params.ratings(:,subI));
%     Generate_params.PriorVar    = var(Generate_params.ratings(:,subI));
    
    % compute prior mean and variance using this_sequence ratings (or prices) instead
    % of whole dataset ratings
%     Generate_params.PriorMean   = mean(Generate_params.seq_vals(sequence,:,subI));
%     Generate_params.PriorVar    = var(Generate_params.seq_vals(sequence,:,subI));
%     Generate_params.PriorMean   = mean(Generate_params.rateseq_vals(sequence,:,subI));
%     Generate_params.PriorVar    = var(Generate_params.rateseq_vals(sequence,:,subI));
% %     
    % ranks for this sequence
    dataList                    = tiedrank(squeeze(Generate_params.rateseq_vals(sequence,:,subI))');
%     dataList                    = tiedrank(squeeze(Generate_params.seq_vals(sequence,:,subI))');
    list.vals                   =  list.allVals;
    
    % here Nick runs the analzeSecretary2021 function. 
    % I am trying to adapt to our best-choice tasks
    % [choiceStop, choiceCont, difVal]  = analyzeSecretaryNick_2021(Generate_params,list);
    
    % Extract params. Hopefully everything about model is now pre-specified
    % (I would update mean and variance directly from Generate_params 
    prior.mu                    = Generate_params.PriorMean + 0;    % prior mean offset is zero unless biased prior model
    prior.sig                   = Generate_params.PriorVar + 0;     % Would a biased variance model be a distinct model?
    if prior.sig < 1; prior.sig = 1; end                            % It can happen randomly that a subject has a low variance and subtracting the bias gives a negative variance. Here, variance is set to minimal possible value.
    prior.kappa                 = Generate_params.kappa;            % prior mean update parameter
    prior.nu                    = Generate_params.nu;
    
    Cs                          = Generate_params.Cs;               % cost-to-sample is set to zero
    
    % bin edges
    minValue                    = Generate_params.binEdges_reward;
    
    list.vals = list.allVals';
    sampleSeries = list.vals;
    N = Generate_params.seq_length;

    [choiceStop, choiceCont, difVal, currentRnk] = computeSecretary(Generate_params, sampleSeries, prior, N, list, Cs,minValue);
    
    num_samples(sequence,subI) = find(difVal<0,1,'first')  % assign output num samples for Bruno model
    
    % assign rank
    ranks(sequence,subI) = dataList(num_samples(sequence,subI));
    % Accumulate action values too so you can compute ll outside this function if needed
    choiceStop_all(sequence, :, subI) = choiceStop;
    choiceCont_all(sequence, :, subI) = choiceCont;
    
end

clear Generate_params choiceStop choiceStop_all choiceCont choiceCont_all difVal currentRnk ranks list N sampleSeries minValue prior dataList
