%% RUN IDEAL OBSERVER %%

% this script tests the ideal observer with the economic task. Here I
% define all the params needed to run the io model. Once that works, I will
% implement it through the main preprocessing script 
% The script operates on one participant at a time. 

% The input/data needed are taken from economic_prepro_v2.m (all data but
% extraction of in-sequences info -- last analysis part) 

% initiate a few parameters
log_or_not      = 0; % log or not?
subjects        = 1; % one subject at a time
subI            = 1; % test subject

%% EXTRACT SUBJECT DATA AND SAVE IN PARAMS %%

% for now we need to extract this_subject phase 1 ratings [400x1], sequences
% [1x40x10] and number of samples [40x1].  
subrates        = allsubs_ratings{1,subI};
subsequences    = allsubs_price_sequences{1,subI};
subsamples      = allsubs_data{1,subI}.samples;
subranks        = allsubs_data{1,subI}.rank;

% 1. get ranks? this seems a bit confusing, because I already have the
% rank of options that each subject chose in allsubs_data.rank struct
% field. Also I don't get the code in line 93 in Nick's code
% for i = 1:length(subsequences)
%     sequenceranks{1,i} = tiedrank(subsequences{1,i})';
% end
% 
% 2 log_or_not 
if log_or_not == 1
    
    Generate_params.ratings(:,sub)          = log(subrates);

    for i = 1:length(subsequences)
        Generate_params.seq_vals(i,:,sub)   = log(subsequences{1,i})';
    end
else
    Generate_params.ratings(:,sub)          = subrates;

    for i = 1:length(subsequences)
        Generate_params.seq_vals(i,:,sub)   = subsequences{1,i}';
    end
end

Generate_params.num_samples(:,sub)          = subsamples;
% how do I save ranks here? Is it the rank of the option that the
% participant chose? 
Generate_params.ranks(:,sub)                = subranks; % ranks may not be needed

% define a few more params 
Generate_params.num_subs                    = size(Generate_params.seq_vals,3);
Generate_params.num_seqs                    = size(Generate_params.seq_vals,1);
Generate_params.seq_length                  = size(Generate_params.seq_vals,2);
Generate_params.num_vals                    = size(Generate_params.ratings,1);
Generate_params.rating_bounds               = [0 100]; % What is min and max of rating scale? (Works for big trust anyway)
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

% model operates one one sequence at a time
for sequence = 1:Generate_params.num_seqs
    
    % squeezing here if more than one subjects
    list.allVals                = squeeze(Generate_params.seq_vals(sequence,:,subI)); % squeezing may not be needed
    Generate_params.PriorMean   = mean(Generate_params.ratings(:,subI));
    Generate_params.PriorVar    = var(Generate_params.ratings(:,subI));
    
    
    
    
end

