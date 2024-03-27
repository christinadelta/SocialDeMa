% % PRE-PROCESSING AND ANALYSIS SCRIPT FOR THE BEADS TASK VERSION 2

% Part of the Optimal Stopping Problems Project

% CREATED: 21/03/2024


%% housekeeping commands

clear all
clc

%% set figure-docking as default 

set(0,'DefaultFigureWindowStyle','docked')

%% INIT LOAD DATA %%
% GET PATHS & DEFINE VARIABLES
% The four next lines (paths) should be changed to your paths 
startpath       = '/Users/christinadelta/githubstuff/rhul_stuff/SocialDeMa/';
datapath        = '/Volumes/DeepSpaceStuff/optimal_stopping_data/data/';
resultspath     = fullfile(datapath, 'beads', 'behav');
croppedpath     = fullfile(startpath, 'analysis', 'beads', 'behav', 'cropped');

task            = 'beads';
subs            = dir(fullfile(resultspath, '*sub*'));
nsubs           = length(subs);
% nsubs           = 5;

totaltrials     = 52; 
conditions      = 2;
condtrials      = totaltrials/conditions;
nmodels         = 2;

% init required variables
avdraws                 = nan(nsubs,1);
avacc                   = nan(nsubs,1);
easy_avdraws            = nan(nsubs,1);
diff_avdraws            = nan(nsubs,1);
easy_avacc              = nan(nsubs,1);
diff_avacc              = nan(nsubs,1);
all_ioacc               = nan(nsubs,2);
all_iodraws             = nan(nsubs,2);

%% PART 1: EXTRACT DATA FROM LOG FILES - RUN IO

for sub = 1:nsubs

    fprintf('loading beads block data\n')  
    subject = subs(sub).name;
    subdir  = fullfile(resultspath,subject);
    fprintf('\t reading data from subject %d\n',sub); 
    
    % extract blocktrial data
    [subsequences,subchoiceVec,all_data, draws_index]    = get_blockdata(subdir,sub,task);
    
    % store in cell for each participant
    allsub_alldata{1,sub}                  = all_data;
    allsub_sequences{1,sub}                = subsequences;
    allsub_choiceVec{1,sub}                = subchoiceVec;
    allsub_drawinfo{1,sub}                 = draws_index; % col1: draw, col2: trial

     %% SPLIT SUB DATA INTO CONDITIONS
    
    for cond = 1:conditions % loop over conditions
        
        tmp                             = find(all_data(:,8) == cond);
        cond_data{sub}{cond}           = all_data((tmp),:);
        clear tmp
        
    end % end of condition loop

    %% AVERAGE PARTICIPANT DRAWS & ACCURACY & CALCULATE POINTS%%
    
    % create a nx1 vector (n=number of participants) with the averaged number
    % of draws for each participant.
    % This vector will be used as a covariate for the individual differences
    % analysis in SPM12.
    
    sub_draws           = all_data(:,5);
    sub_acc             = all_data(:,7);
    avdraws(sub,1)     = mean(sub_draws);
    avacc(sub,1)       = mean(sub_acc);
    
    clear sub_draws sub_acc
    
    % average draws and acc for each condition
    for cond = 1:conditions
        
        tmp_cond                = cond_data{1,sub}{1,cond};
        cond_draws              = tmp_cond(:,5);
        cond_acc                = tmp_cond(:,7);
        
        if cond == 1
            easy_avdraws(sub,1) = mean(cond_draws);
            easy_avacc(sub,1)   = mean(cond_acc);
        else
            diff_avdraws(sub,1) = mean(cond_draws);
            diff_avacc(sub,1)   = mean(cond_acc);
            
        end
        
        % we will now calculate participant points. we can use this to se how each participant performed and to compare
        % each participant with their corresponding model instant 
    
        % we will need: cost_correct, cost_wrong, cost_to_sample, numdraws, acc 
        allsub_points(sub, cond) = (sum(cond_acc==1)*10) + (sum(cond_acc==0)*-10) + (sum(cond_draws)*-0.25);
        
    end 

    %% RUN IDEAL OBSERVER 
    
    % define parameters of the ideal observer
    R.alpha             = 1;            % softmax stochasticity parameter (for fitting to human behaviour) - this is not needed here
    R.error             = -10;          % cost for being wrong
    R.correct           = 10;           % reward for being correct
    R.q                 = [0.8 0.6];    % proportion of the majority value in sequence (60:40 split in this case)
    R.sample            = -0.25;        % the cost to sample

    % loo over conditions
    for cond = 1:conditions

        thiscond_data                   = cond_data{1,sub}{1,cond};
        thiscond_seq                    = subsequences{1,cond};
        R.cond                          = cond;
        io_output                       = run_POMDP_BeadsIo(R,thiscond_seq,thiscond_data);

        % extract output
        allsub_io_output{1,sub}{1,cond} = io_output;
        all_ioacc(sub,cond)             = io_output.accuracy;
        all_iodraws(sub,cond)           = io_output.draws;  

        clear io_output
    end 

end % end of subjects loop

%% PART 2: RUN STATISTICS (all statistcally significant)

% add draws and performance in one vec
all_acc             = [easy_avacc diff_avacc];
all_draws           = [easy_avdraws diff_avdraws];

% make a struct with all the vectors needed for the analysis
anova_struct        = struct('all_draws', all_draws, 'all_acc', all_acc,...
    'all_ioacc', all_ioacc, 'all_iodraws', all_iodraws);

% run anovas and pairwise comparisons 
output_struct_one   = runBehav_stats(nsubs, anova_struct); % output will be used for plotting 

clear anova_struct

%% PART 3: FIT BETA+CS MODEL TO PARTICIPANT DATA 

% define fixed parameters used in all models 
R.error             = -10;          % cost for being wrong
R.correct           = 10;           % reward for being correct
R.q                 = [0.8 0.6];    % proportion of the majority value in sequence (60:40 split in this case)
R.difference        = -20;
R.Cs                = -0.25;
R.beta              = 3;

for sub = 1:nsubs

    % define free parameters
    R.initsample        = R.Cs;
    R.initbeta          = R.beta;

    subdata             = cond_data{1,sub};         % all data matrix
    sub_choicesvec      = allsub_choiceVec{1,sub};  % subject choices (two 26 by 3 matricies)
    sub_sequence        = allsub_sequences{1,sub};  % sequences that this-subject was presented with




end% end of subject loop

