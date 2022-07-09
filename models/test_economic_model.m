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
analyse_value_positions = 1;















