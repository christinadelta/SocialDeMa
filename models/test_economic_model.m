%%% Test model fitting and ideal observer for Best-Choice task

simulate_stimuli        = 1; % If 1, randomly generates phase 1 distributions and sequences from them and saves info in Generate_params
check_params            = 0; % fit the same model that created the data and output estimated parameters
make_est_model_data     = 0;
use_file_for_plots      = 1;
make_plots              = 1; % if 1, plots the results
all_draws_set           = 1; % You can toggle how the ll is computed here for all models at once if you want or go bvelow and set different values for different models manually in structure
log_or_not              = 0; % I'm changing things so all simulated data is logged at point of simulation (==1) or not

% models:
% 1: cutoff
% 2: Cs
% 3: dummy
% 4: BV
% 5: BR
% 6: BPM
% 7: Opt
% 8: BPV

% These are now what v2 called model identifiers
do_models               = [1 2 4 5 7 ]; % Applies at the moment to both make_model_data and check_params;
comment                 = sprintf('facesonlineSubsLog%d',log_or_not); %The filename will already fill in basic parameters so only use special info for this.
%These correspond to identifiers (not configured implementations like in v2) in the v3_sweep version
model_names             = {'Cut off' 'Cs' 'IO' 'BV' 'BR' 'BPM' 'Opt' 'BPV' }; % IO is a placeholder, don't implement
num_model_identifiers   = size(model_names,2);
subjects                = 1;
IC                      = 2; % either AIC (1) or BIC (2)
analyze_value_positions = 1;
outpath                 = '/Users/christinadelta/githubstuff/rhul_stuff/SocialDeMa/models';
counter                 = 0; % this is used in Nick's code, in the subject's loop (may not use it, but will leave it for now)

% loop over subjects
for sub = 1: subjects
    
    % for now we need to extract this_subject phase 1 ratings [400x1], sequences
    % [1x40x10] and number of samples [40x1]. Here 1=number of subjects 
    subrates        = allsubs_ratings{1,sub};
    subsequences    = allsubs_sequences{1,sub};
    subsamples      = allsubs_data{1,sub}.samples;
    subranks        = allsubs_data{1,sub}.rank;

    % 1. get ranks? this seems a bit confusing, because I already have the
    % rank of options that each subject chose in allsubs_data.rank struct
    % field. Also I don't get the code in line 93 in Nick's code
%     for i = 1:length(subsequences)
%         sequenceranks{1,i} = tiedrank(subsequences{1,i})';
%     end
    
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
    Generate_params.ranks(:,sub)                = subranks;
    
    % load the params struct with all the info needed
    Generate_params.analyze_value_positions     = analyze_value_positions;  % make psychometric plots if set to 1
    Generate_params.IC                          = IC; % AIC (1) or BIC (2) correction?
    Generate_params.log_or_not                  = log_or_not;
    Generate_params.all_draws_set               = all_draws_set; %% what is that exactly?
    Generate_params.do_models_identifiers       = do_models;
    Generate_params.num_subs                    = size(Generate_params.seq_vals,3);
    Generate_params.num_seqs                    = size(Generate_params.seq_vals,1);
    Generate_params.seq_length                  = size(Generate_params.seq_vals,2);
    Generate_params.num_vals                    = size(Generate_params.ratings,1);
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
    
    
    %% SET UP MODELS %%
    
    %We build off of (Or rather use part of) the infrastructure I created for configuring models in
    %Param_recover*.m. It involves a template list of default parameters
    %values that is then repmatted into a parameter*model type matrix of
    %default parameters. Then a matching parameter*model type free_parameters
    %matrix marks which parameters to fit down below.These matrices are
    %then used to populate separate model fields in Generate_params, which
    %is then dropped into the estimation function.
    
    %Make the template parameter list
    opt_rule                    = ceil(exp(-1) * Generate_params.seq_length);  % 37% rule cutoff
    
    model_template.identifier   = 2;                % row 1 in param_config, 1:CO 2:IO 3:Cs 4:BV 5:BR 6:BP 7:optimism 8:BPV
    model_template.kappa        = 2;                % row 2 in param_config
    model_template.nu           = 1;                % row 3
    model_template.cutoff       = opt_rule;         % row 4, initialised to optimal 37%
    model_template.Cs           = 0;                % row 5, intiialised to optimal no cost to sample
    model_template.BVslope      = 0.2;              % row 6, intialised to 1 (like a threshold)
    model_template.BVmid        = 55;               % row 7, initialised to halfway through the rating scale (can't be used with log)
    model_template.BRslope      = 1;                % row 8
    model_template.BRmid        = 55;               % row 9
    model_template.BP           = 0;                % row 10
    model_template.optimism     = 0;                % row 11
    model_template.BPV          = 0;                % row 12
    model_template.log_or_not   = log_or_not;       % 1 = log transform (normalise) ratings  %row 13 (This is already legacy - log or not is now controlled by switch at start of programme and simulated data was logged before reaching this point).
    model_template.all_draws    = all_draws_set;    % 1 = use all trials when computing ll instead of last two per sequence.   %row 14
    model_template.beta         = 1;                % Just for parameter estimation.
    model_template.name         = 'template';
    
    % Correct the starting parameters that are in units of ratings
    if log_or_not == 1
        model_template.BVmid    = log(model_template.BVmid);
        model_template.BRmid    = log(model_template.BRmid);
        
        % BP & BVP would be in log units too but can't take negative or
        % 0th values so it's best to manually set their starting params
        % and let the estimation find the best value for the context.
        % That means fix them manually if you want them to be logged
    end
    
    % Repmat the template to create a column for each model. For now, we are
    % doing all possible models, not the ones specified in do_models. We'll
    % reduce this matrix to just those below.
    identifiers                 = 1:num_model_identifiers;
    num_cols                    = num_model_identifiers;
    param_config_default        = [ ...
        identifiers;                                    % row 1: identifiers,  1:CO 2:IO 3:Cs 4:BV 5:BR 6:BP 7:optimism 8:BPV
        repmat(model_template.kappa,1,num_cols);        % row 2: kappa
        repmat(model_template.nu,1,num_cols);           % row 3: nu
        repmat(model_template.cutoff,1,num_cols)        % row 4: cutoff
        repmat(model_template.Cs,1,num_cols);           % row 5: Cs
        repmat(model_template.BVslope,1,num_cols);      % row 6: BV slope
        repmat(model_template.BVmid,1,num_cols);        % row 7: BV mid
        repmat(model_template.BRslope,1,num_cols);      % row 8: BR slope
        repmat(model_template.BRmid,1,num_cols);        % row 9: BR mid
        repmat(model_template.BP,1,num_cols);           % row 10: prior mean offset (BP)
        repmat(model_template.optimism,1,num_cols);     % row 11: optimism
        repmat(model_template.BPV,1,num_cols);          % row 12: prior variance offset (BPV)
        % repmat(model_template.log_or_not,1,num_cols); % row 13: log or not (at the moment not to be trusted)
        repmat(model_template.all_draws,1,num_cols);    % row 14: all draws
        repmat(model_template.beta,1,num_cols);         % row 15: beta
        ];
    
    % Mark which are free/to be estimated
    free_parameters         = zeros(size(param_config_default));
    free_parameters(4,1)    = 1;    % Model indicator 1, parameter 4: Cut off
    free_parameters(5,2)    = 1;    % Model indicator 2, parameter 5: Cs
    free_parameters(7,4)    = 1;    % Model indicator 4, parameter 7: BV
    free_parameters(9,5)    = 1;    % Model indicator 5, parameter 9: BR
    free_parameters(10,6)   = 1;    % Model indicator 6, parameter 10: BPM
    free_parameters(11,7)   = 1;    % Model indicator 7, parameter 11: Opt
    free_parameters(12,8)   = 1;    % Model indicator 8, parameter 12: BPV
    
    % Now reduce matrices to just those in do_models
    % In Param_recover*.m we had distrinctions between model instantiations
    % and
    param_config_default                    = param_config_default(:,do_models);
    free_parameters                         = free_parameters(:,do_models);
    
    % Save your work into struct
    Generate_params.num_models              = numel(do_models);
    Generate_params.param_config_default    = param_config_default;
    Generate_params.free_parameters_matrix  = free_parameters;
    Generate_params.comment                 = comment;
    Generate_params.outpath                 = outpath;
    analysis_name                           = sprintf(...
        'out_new_ll%d_%s_'...
        , Generate_params.all_draws_set ...
        , Generate_params.comment ...
        );
    Generate_params.analysis_name           = analysis_name;
    outname                                 = [analysis_name char(datetime('now','format','yyyyddMM')) '.mat'];
    Generate_params.outname                 = outname;
    
    







end % 
