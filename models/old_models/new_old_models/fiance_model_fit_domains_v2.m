
function [] = fiance_model_fit_domains_v2;

%Don't forget, if you want to process cluster data, the individual beta
%files are in C:\matlab_files\fiance\online_domains_01\fitted_datafiles and
%you use domains_process_cluster_data.m to make new files with the best
%beta from the grid search. Then you need to use this function with
%run_fitted data on the parameters from the best beta first before plotting.

%v2 just starts to save result files in new folder so they quit
%accumulating in my main directory. I still need the old version for now
%though so am saving it.

estimate_param = 0; %If 1 it will estimate parameters of the models specified in param_to_fit. If 0, it will open pre-saved files for the specified models and continue analyses of those (e.g., draws/ranks, comparisons with hum subjects)
param_to_fit = [1:3 5];  %1==Cs; 2==biased values (slope,mid). 3==ideal observer(no params but fitted beta) 4 == (no longer supported) 5 == biased prior 6 = cut off; Can use vector of model numbers and it'll run them all
run_fitted_models = 0;  %If 1, will calculate results of models using its previously fitted parameters
make_plots = 1; %If 1, will open files containing the results of the fitted models and plot them.
analyze_value_positions = 1;    %analyses and makes sequential threshold plots of fitted data (This won't work unless make_plots also equals 1)
experiment = 1; %1: online holidays 2: online foods 3: online  open day
IC = 1; %1 if AIC, 2 if BIC
two_params = 1;
betas_to_test = [1:20];
log_or_not = 1; %1 log transforms values before decision making

no_params = [2 3 1 1 2 2]; %Used for IC correction
nbins_psi = 6;  %for serial position analysis
nbins_reward = 6;  %For assigning reward value to outcomes in model

num_models = numel(param_to_fit);

if log_or_not == 1;
    %from cluster,uses log and alldraws on indiv. subs
            model_strs = {'Cs_lall' 'BV_lall' 'io_lall' 'X' 'BP_lall' 'CO_lall'};    %
%     model_strs = {'Cs_lnss' 'BV_lss' 'io_lss' 'X' 'BP_lss' 'CO_lss'}; 
else
        model_strs = {'Cs_ss' 'BV_ss' 'io_ss' 'X' 'BP_ss' 'CO_ss'};    %
%     model_strs = {'Cs_n' 'BV_n' 'io_grid' 'BVo_n' 'BP_n' 'CO_n'};
end;

for model = 1:num_models;
    
    %For fitting, I'm going to run and save each model separately without
    %accumulating anything over models. This is backwards compatible with
    %previous versions and makes model analysis more robust to differences
    %in the order/cuircumstances underwhich models were fitted.
    clear mparams lla list sub_choice_trial sub_choice_rank dataPrior difVal;
    
    %now returns (preserving legacy variable names):
    %"mean ratings" which is actually 90*num_subs lists of phase 1 ratings
    %seq_vals, which is 6*8*num_subs lists of sequence values and
    %output, which is now 6*num_subs number of subject draws for each sequence
    [mean_ratings seq_vals output] = get_sub_data(experiment);
    
    %Now things become tangled because the model fitting needs the subject
    %loop to be inside of f_fitparams.m so that beta can be fit outside it
    
    if estimate_param == 1;
        
        param_initial_beta = 7;
        
        %just jump straight into outer, beta fit
        disp(sprintf('fitting beta model %d',param_to_fit(model)));
        
        %fit_outer and fminsearch will fit beta to group and other free
        %parameters to individual subjects but for now it only works for
        %ideal observer (which has no individual participant free
        %paraneters) because the free parameters that are computed and
        %still trapped inside fminsearch and I never figured out how to
        %return them to this base function to be saved.
        %         options = optimset('Display','iter');
        %         [param_fitted_beta ll_beta exitflag(model) search_out(model)] = ...
        %             fminsearch(@(param_initial_beta) fit_outer(param_initial_beta, param_to_fit(model), mean_ratings, seq_vals, output, nbins_reward, experiment), param_initial_beta, options, log_or_not);
        %
        %fit_outer grid searches integer values of beta for the one whos
        %ll's summed over individual model fits is minimised and iot
        %returns the fitted individual participant parameters
        [param_fitted_beta ll_beta mparams mlla] = fit_outer_grid(param_to_fit(model), mean_ratings, seq_vals, output, nbins_reward, experiment, betas_to_test, log_or_not);
        
        disp(['saving C:\matlab_files\fiance\online_domains_01\fitted_datafiles\fit_domains_' model_strs{param_to_fit(model)} sprintf('_params_exp%d.mat',experiment)]);
        save(['C:\matlab_files\fiance\online_domains_01\fitted_datafiles\fit_domains_' model_strs{param_to_fit(model)} sprintf('_params_exp%d.mat',experiment)]);
        
    end; %estimate params?
    
    if run_fitted_models == 1;
        
        %returns subject*sequences matrices of numbers of draws and ranks
        [fitted_choice_trial fitted_choice_rank] = get_a_models_performance(experiment, param_to_fit(model), model_strs, mean_ratings, seq_vals, nbins_reward, log_or_not);
        
        disp(['saving C:\matlab_files\fiance\online_domains_01\fitted_datafiles\fitted_domains_' model_strs{param_to_fit(model)} sprintf('_params_exp%d.mat',experiment)]);
        save(['C:\matlab_files\fiance\online_domains_01\fitted_datafiles\fitted_domains_' model_strs{param_to_fit(model)} sprintf('_params_exp%d.mat',experiment)]);
        
    end; %run fitted models?
    
end;    %loop through models

if make_plots == 1;
    
    plot_data(experiment, model_strs, param_to_fit, mean_ratings, output, seq_vals, IC, two_params, no_params, analyze_value_positions, nbins_psi, log_or_not )
    
    
end;    %make plots?

disp('audi5000');

function  [param_fitted_beta ll_beta mparams mlla] = fit_outer_grid(which_model, all_ratings, all_seq_vals, all_output, nbins_reward, experiment, betas_to_test, log_or_not)

ll_beta = 0;

num_subs = size(all_ratings,2);
max_subs = num_subs;

for beta_val = 1:numel(betas_to_test);
    
    param_initial_beta = betas_to_test(beta_val);
    
    for num_subs_found = 1:max_subs;
        
        
        %peel off this sub
        seq_vals = all_seq_vals(:,:,num_subs_found);
        mean_ratings = all_ratings(:,num_subs_found);   %returns to original naming convention, at risk of being confusing
        
        %set this sub's prior
        if log_or_not == 1;
            binEdges_reward = linspace(min(log(mean_ratings(:))),max(log(mean_ratings(:))),nbins_reward+1);   %organise bins by min and max
            dataPrior.mean = mean(log(mean_ratings));   %take the log to normalize
            dataPrior.var = var(log(mean_ratings));
            dataPrior.range = [min(min(log(seq_vals))) max(max(log(seq_vals)))];   %used for BV normilisation
        else
            binEdges_reward = linspace(min(mean_ratings(:)),max(mean_ratings(:)),nbins_reward+1);   %organise bins by min and max
            dataPrior.mean = mean(mean_ratings);  
            dataPrior.var = var(mean_ratings);
            dataPrior.range = [min(min(seq_vals)) max(max(seq_vals))];   %used for BV normilisation
        end;
        distOptions = 0;
        dataPrior.kappa = 2;
        dataPrior.nu = 1;
        %     options = optimset('Display','iter','MaxFunEvals', 5000, 'TolFun', 0.001);
        
        %now loop through the sequences, just to set up the data to be passed
        %into fit_inner
        for sequence=1:size(seq_vals,1);
            
            %ranks for this sequence
            dataList = tiedrank(seq_vals(sequence,:)');
            
            %Subject draws and ranks
            sub_choice_trial(num_subs_found,sequence) = all_output(sequence,num_subs_found);  %just transferring draws from new variable name back to old one
            sub_choice_rank(num_subs_found,sequence) = dataList(sub_choice_trial(num_subs_found,sequence)); %ranks
            
            %format data for model fitting
            %initialze list
            if log_or_not == 1;
                list(sequence).allVals = log(seq_vals(sequence,:));
                list(sequence).vals = log(seq_vals(sequence,1:sub_choice_trial(num_subs_found,sequence)));
            else;
                list(sequence).allVals = seq_vals(sequence,:);
                list(sequence).vals = seq_vals(sequence,1:sub_choice_trial(num_subs_found,sequence));
            end;
            list(sequence).optimize = 0;
            list(sequence).flip = 1;
            list(sequence).length = size( list(sequence).allVals, 2 );
            
        end;    %loop through sequences (will start new sequence loop inside fit_inner to run model on each sequence, this loop here just compiles the data to be passed to fit_inner
        
        if which_model == 1 %cost to sample fit
            
            params = [-0.03];   %initialise Cs & beta
            
        elseif which_model == 2;
            
%             %initialise parameter search using psychometric function values
%             if experiment == 1; %av
%                 load('C:\matlab_files\fiance\power_param_study2_6bins_indsub1.mat','parameters');  %pre-fit and used in publication
%             elseif experiment == 2; %matchmaker
%                 load('C:\matlab_files\fiance\power_param_study9_6bins_indsub1.mat','parameters');  %pre-fit and used in publication
%             elseif experiment == 3; %trustworthiness
%                 load('C:\matlab_files\fiance\power_param_study13_6bins_indsub1.mat','parameters');  %pre-fit and used in publication
%             end;
%             
%             params = [parameters(:,num_subs_found)']; %initialise with this subs' data and beta=5;

            params = [1 50];
            
        elseif which_model == 3 | which_model == 4;    %ideal observer and biased values pre-fit will fit only beta (prefit logistic paraneters fixed inside sectretary progtramme)
            
            params = [0]; %I'll just pass in a dummy parameter which does affect computations any so I can use the same function for testing.
            
        elseif which_model == 5;    %biased prior
            
            params = [0.5]; %initialise difference from mean and beta too
            
        elseif which_model == 6;    %cut-off
            
            params = (1/exp(1))*list(1).length; %initialise cut off to 37% rule
            
        end;
        
        %     disp(sprintf('fitting model %d to subject %d', which_model, num_subs_found));
        if num_subs_found == 1;
            fprintf('beta %d: subject 1, ',param_initial_beta);
        else
            fprintf('%d, ',num_subs_found);
        end;
        
        %     options = optimset('Display','iter');
        [temp_params(num_subs_found,:), lla(num_subs_found,beta_val), exitflag(num_subs_found), search_out(num_subs_found)] = ...
            fminsearch(@(params) fit_inner(params, dataPrior, list, distOptions,binEdges_reward,which_model,experiment,num_subs_found,param_initial_beta), params);
        
    end;    %loop through subjects
    
    temp_params_beta(:,:,beta_val) = temp_params;
    temp_ll_beta(beta_val) = sum(lla(:,beta_val));
    
    fprintf('... done,  ll: %3.2f\n', temp_ll_beta(beta_val));
    
end;    %loop through betas

[ll_beta I] = min(temp_ll_beta); %find the summed log likelihood of beta with the lowest sum log likelihood out of all of the betas tried
param_fitted_beta = betas_to_test(I);   %take out the beta value that corresponds to the lowest summed log likelihood
mparams = squeeze(temp_params_beta(:,:,I)); %take out the individual subject parameters that have the lowest summed log likelihoods
mlla = lla(:,I);    %take out the individual subject log likelihoods that have the lowest sum out of all betas tried


%%

function plot_data(experiment, model_strs, param_to_fit, mean_ratings, output, seq_vals, IC, two_params, no_params, analyze_value_positions, nbins_psi , log_or_not);

%Ok I know I just closed the model loop but now, to preserve modularity, I
%will immediateky re-open a new model loop, so I can make plots
%independently of fitting the models or running the fitted models.

%variable renaming to new scheme and remake sub_choice_rank
sub_choice_trial = output';
for sub=1:size(seq_vals,3);
    
    if log_or_not == 1;
        binEdges_psi(sub,:) = linspace(min(log(mean_ratings(:,sub))),max(log(mean_ratings(:,sub))),nbins_psi+1);
    else;
        binEdges_psi(sub,:) = linspace(min(mean_ratings(:,sub)),max(mean_ratings(:,sub)),nbins_psi+1);
    end;
    
    for seq=1:size(seq_vals,1);
        dataList = tiedrank(seq_vals(seq,:,sub));
        sub_choice_rank(sub,seq) = dataList(sub_choice_trial(sub,seq));
    end;    %seq
end;    %subs

%regretting the modular functions now
num_models = numel(param_to_fit);

plot_cmap = hsv(size(model_strs,2)+1);  %models + subjects
f_a = 0.1; %face alpha
sw = 0.5;  %ppoint spread width
font_size = 12;

%first make plot of draws
h1 = figure; set(gcf,'Color',[1 1 1]);

%Before embarking on models, plot participant draws in first bar position
subplot(2,2,1); hold on;

%sub datapoints
handles = plotSpread(mean(sub_choice_trial,2), ...
    'xValues',0,'distributionColors',plot_cmap(1,:),'distributionMarkers','.', 'spreadWidth', sw);

bar(0,mean(mean(sub_choice_trial,2)), ...
    'FaceColor',plot_cmap(1,:),'FaceAlpha',f_a,'EdgeColor',[0 0 0] );

ylim([0 8]);
ylabel('Samples to decision');
set(gca,'XTick',[1:num_models+1],...
    'YTick',[2:2:8], ...
    'FontSize',font_size,'FontName','Arial');
box off;


%Before embarking on models, plot participant ranks in first bar position
subplot(2,2,2); hold on;

%sub datapoints
handles = plotSpread(mean(sub_choice_rank,2), ...
    'xValues',0,'distributionColors',plot_cmap(1,:),'distributionMarkers','.', 'spreadWidth', sw);

bar(0,mean(mean(sub_choice_rank,2)), ...
    'FaceColor',plot_cmap(1,:),'FaceAlpha',f_a,'EdgeColor',[0 0 0] );

ylim([0 8]);
ylabel('Samples to decision');
set(gca,'XTick',[1:num_models+1],...
    'YTick',[2:2:8], ...
    'FontSize',font_size,'FontName','Arial');
box off;

%     %Accumulate choice data
%     if experiment == 3; %something is wrong with subject 18;
%         all_choice_trial(:,:,1) = sub_choice_trial([1:17 19 20],:);
%     else
all_choice_trial(:,:,1) = sub_choice_trial;
%     end;
%
for model = 1:num_models;
    
    fprintf(' ');
    
    clear fitted_choice_trial fitted_choice_rank;
    load(['C:\matlab_files\fiance\online_domains_01\fitted_datafiles\fit_domains_' model_strs{param_to_fit(model)} sprintf('_params_exp%d.mat',experiment)],'mlla');
    load(['C:\matlab_files\fiance\online_domains_01\fitted_datafiles\fitted_domains_' model_strs{param_to_fit(model)} sprintf('_params_exp%d.mat',experiment)],'fitted_choice_trial', 'fitted_choice_rank');
    
    lla = mlla';
    all_choice_trial(:,:,model+1) = fitted_choice_trial;
    
    %Model samples
    figure(h1); subplot(2,2,1);
    %sub datapoints
    handles = plotSpread(nanmean(fitted_choice_trial,2), ...
        'xValues',model,'distributionColors',plot_cmap(param_to_fit(model)+1,:),'distributionMarkers','.', 'spreadWidth', sw);
    
    bar(model,nanmean(nanmean(fitted_choice_trial,2)), ...
        'FaceColor',plot_cmap(param_to_fit(model)+1,:),'FaceAlpha',f_a,'EdgeColor',[0 0 0] );
    
    set(gca,'XTick',[1:num_models],'fontSize',font_size,'FontName','Arial',...
        'XLim',[-1 num_models+2],'YLim',[0 8]);
    ylabel('Samples to decision');
    xlabel('Model');
    
    %model ranks
    figure(h1); subplot(2,2,2);
    %sub datapoints
    handles = plotSpread(mean(fitted_choice_rank,2), ...
        'xValues',model,'distributionColors',plot_cmap(param_to_fit(model)+1,:),'distributionMarkers','.', 'spreadWidth', sw);
    
    bar(model,nanmean(nanmean(fitted_choice_rank,2)), ...
        'FaceColor',plot_cmap(param_to_fit(model)+1,:),'FaceAlpha',f_a,'EdgeColor',[0 0 0] );
    
    set(gca,'XTick',[1:num_models],'fontSize',font_size,'FontName','Arial',...
        'XLim',[-1 num_models+2],'YLim',[0 8]);
    ylabel('Rank of chosen option');
    xlabel('Model');
    
    %model likelihoods
    figure(h1); subplot(2,4,5); hold on;
    
    %sub datapoints
    handles = plotSpread(lla', ...
        'xValues',model,'distributionColors',plot_cmap(param_to_fit(model)+1,:),'distributionMarkers','.', 'spreadWidth', sw);
    
    bar(model,nanmean(lla), ...
        'FaceColor',plot_cmap(param_to_fit(model)+1,:),'FaceAlpha',f_a,'EdgeColor',[0 0 0] );
    
    set(gca,'XTick',[1:num_models],'fontSize',font_size,'FontName','Arial',...
        'XLim',[0 num_models+2]);
    %,'YLim',[15 40]);
    ylabel('log-likelihood');
    xlabel('Model');
    
    %Model IC
    
    if IC == 1; %If AIC (per participant)
        IC_pps = 2*no_params(param_to_fit(model)) + 2*lla;
        IC_sum = nansum(IC_pps);
        %             IC_sum = 2*no_params(param_to_fit(model)) + 2*nansum(lla);
        a_label = 'AIC';
        IC_ylims = [800 1350];
    elseif IC == 2; %If BIC (per participant)
        IC_pps = no_params(param_to_fit(model))*log(28) + 2*lla;
        IC_sum = nansum(IC_pps);
        %             IC_sum = no_params(param_to_fit(model))*log(numel(lla)*28) + 2*nansum(lla);
        a_label = 'BIC';
        IC_ylims = [750 1250];
    end;
    
    %Plot of IC averages
    subplot(2,4,6); hold on;
    
    handles = plotSpread(IC_pps', ...
        'xValues',model,'distributionColors',plot_cmap(param_to_fit(model)+1,:),'distributionMarkers','.', 'spreadWidth', sw);
    
    bar(model,nanmean(IC_pps), ...
        'FaceColor',plot_cmap(param_to_fit(model)+1,:),'FaceAlpha',f_a,'EdgeColor',[0 0 0] );
    
    set(gca,'XTick',[1:num_models],'fontSize',font_size,'FontName','Arial',...
        'XLim',[0 num_models+2]);
    ylabel(a_label);
    xlabel('Model');
    
    IC_pps_all_models(:,model) = IC_pps';
    
    %sum of of IC measure
    subplot(2,4,7); hold on;
    
    bar(model,IC_sum, ...
        'FaceColor',plot_cmap(param_to_fit(model)+1,:),'FaceAlpha',f_a,'EdgeColor',[0 0 0] );
    ylabel(sprintf('%s sum',a_label));
    xlabel('Model');
    
    set(gca,'XTick',[1:num_models],'fontSize',font_size,'FontName','Arial',...
        'XLim',[0 num_models+2]);
    %,'YLim',[IC_ylims(1) IC_ylims(2)]);
    ylabel(a_label);
    xlabel('Model');
    
    IC_sum_all_models(:,model) = IC_pps';
    
end;    %loop through models

%Now that model loop is over and data on all models is accumulated, do
%analyses/plots comparing models ...

%ttests on model averages
if num_models > 1;
    
    %run and plot ttests on IC averages
    pairs = nchoosek(1:num_models,2);
    num_pairs = size(pairs,1);
    [a In] = sort(diff(pairs')','descend');  %lengths of connecting lines
    line_pair_order = pairs(In,:);    %move longest connections to top
    
    %Where to put top line?
    y_inc = 2;
    ystart = max(max(IC_pps_all_models)) + y_inc*num_pairs;
    line_y_values = ystart:-y_inc:0;
    
    for pair = 1:num_pairs;
        
        %run ttest this pair
        [h IC_pp_pvals(pair) ci stats] = ttest(IC_pps_all_models(:,line_pair_order(pair,1)), IC_pps_all_models(:,line_pair_order(pair,2)));
        
        %plot result
        subplot(2,4,6); hold on;
        set(gca,'Ylim',[0 ystart]);
        
        if IC_pp_pvals(pair) < 0.05/size(pairs,1);  %multiple comparison corrected
            
            plot([line_pair_order(pair,1) line_pair_order(pair,2)],...
                [line_y_values(pair) line_y_values(pair)],'LineWidth',2,'Color',[0 0 0]);
            
        end;    %Do line on plot?;
        
    end;    %loop through ttest pairs
    
end;    %do I have more than 1 models, for ttests?


%winning models
subplot(2,4,8); hold on;

%winning models
[a pps_indices] = min(IC_sum_all_models');

for model = 1:num_models;
    
    bar(model,numel(find(pps_indices==model)), ...
        'FaceColor',plot_cmap(param_to_fit(model)+1,:),'FaceAlpha',f_a,'EdgeColor',[0 0 0] );
    
end;    %models

set(gca,'XTick',[1:num_models],'fontSize',font_size,'FontName','Arial',...
    'XLim',[0 num_models+2]);
ylabel('Frequency');
xlabel('Model');

if analyze_value_positions == 1;
    
    if log_or_not == 1;
        analyze_value_position_functions(log(seq_vals),all_choice_trial,plot_cmap,binEdges_psi,model_strs,param_to_fit,two_params);
    else
        analyze_value_position_functions(seq_vals,all_choice_trial,plot_cmap,binEdges_psi,model_strs,param_to_fit,two_params);
    end;
    
end;    %make threshold by serial position plot







function [fitted_choice_trial fitted_choice_rank] = get_a_models_performance(experiment, which_model, model_strs, all_ratings, all_seq_vals, nbins_reward,log_or_not);

%returns subject*sequences matrices of numbers of draws and ranks

%note: What is called which_model here is param_to_fit(model) outside this programme,
%at the base level

%get fitted parameters for this model
load(['C:\matlab_files\fiance\online_domains_01\fitted_datafiles\fit_domains_' model_strs{which_model} sprintf('_params_exp%d.mat',experiment)],'mparams');

num_subs = size(all_ratings,2);

for num_subs_found = 1:num_subs;
    
    disp(sprintf('Running already fitted model %d for subject %d to get performance data', which_model, num_subs_found));
    
    try
        use_params = mparams(num_subs_found,:);
    catch
        use_params = [];   %if no parameters get loaded (e.g., ideal observer) this blank will descend into analyzeSecretary
    end;
    
    %set up this subject's model
    %peel off this sub
    seq_vals = all_seq_vals(:,:,num_subs_found);
    mean_ratings = all_ratings(:,num_subs_found);   %returns to original naming convention, at risk of being confusing
    
    %Need for assigning reward
    if log_or_not == 1;
        binEdges_reward = linspace(min(log(mean_ratings(:))),max(log(mean_ratings(:))),nbins_reward+1);   %organise bins by min and max
    else
        binEdges_reward = linspace(min(mean_ratings(:)),max(mean_ratings(:)),nbins_reward+1);   %organise bins by min and max
    end;
    
    %set this sub's prior
    distOptions = 0;
    if log_or_not == 1;
        dataPrior.mean = mean(log(mean_ratings));   %take the log to normalize
        dataPrior.var = var(log(mean_ratings));
        dataPrior.range = [min(min(log(seq_vals))) max(max(log(seq_vals)))];   %used for BV normilisation
    else;
        dataPrior.mean = mean(mean_ratings);   
        dataPrior.var = var(mean_ratings);
        dataPrior.range = [min(min(seq_vals)) max(max(seq_vals))];   %used for BV normilisation
    end;
    dataPrior.kappa = 2;
    dataPrior.nu = 1;
    %     options = optimset('Display','iter','MaxFunEvals', 5000, 'TolFun', 0.001);
    
    for sequence = 1:size(seq_vals,1);
        
        %ranks for this sequence
        dataList = tiedrank(seq_vals(sequence,:)');
        
        %format data for model fitting
        %initialze list
        if log_or_not == 1;
            list(sequence).allVals = log(seq_vals(sequence,:));
        else;
            list(sequence).allVals = seq_vals(sequence,:);
        end;
        
        list(sequence).optimize = 0;
        list(sequence).flip = 1;
        list(sequence).vals =  list(sequence).allVals;
        list(sequence).length = size( list(sequence).allVals, 2 );
        
        
        if which_model == 6;    %cut off rule;
            
            %put ones for every candidate choice
            this_seq_vals = list(sequence).allVals;
            choiceStop = zeros(size(this_seq_vals));
            choiceStop(1,(find( this_seq_vals(ceil(mparams(1))+1:end) > max(this_seq_vals(1:ceil(mparams(1))))))+ceil(mparams(1))) = 1;
            choiceStop(1,list(sequence).length) = 1;
            %put zeros for every non-candidate
            choiceCont = ones(size(this_seq_vals));
            choiceCont(1,find(choiceStop==1))=0;
            difVal(sequence,:,num_subs_found) = choiceCont-choiceStop;
            
        else;   %anything else
            [choiceStop(sequence,:,num_subs_found), choiceCont(sequence,:,num_subs_found), difVal(sequence,:,num_subs_found), currentRnk(sequence,:,num_subs_found), winnings(sequence,num_subs_found)]  = ...
                analyzeSecretaryNick_2020(dataPrior, list(sequence), use_params, distOptions,binEdges_reward,which_model,experiment,num_subs_found);
            %                     analyzeSecretary_brunoversiontest_2020(dataPrior, list(sequence), mparams(num_subs_found,:), distOptions,param_to_fit(model));
        end;
        
        temp = find(difVal(sequence,:,num_subs_found)<0);  %which choices were better to accept than reject?
        if isempty(temp);
            disp(sprintf('problem computing model behaviour for model %d num_subs_found: %d sequence: %d',num_subs_found,sequence));
            fitted_choice_trial(num_subs_found,sequence) = NaN;
            fitted_choice_rank(num_subs_found,sequence) = NaN;
        else
            fitted_choice_trial(num_subs_found,sequence) = temp(1);  %on what trial was the first time value of smpling dropped below value of take?
            dataList = tiedrank(list(sequence).allVals); %get the ranks for this sequence
            fitted_choice_rank(num_subs_found,sequence) = dataList(fitted_choice_trial(num_subs_found,sequence));
        end;
        
    end;    %loop through sequences
    
    
end;    %subject loop






function  ll_beta = fit_outer(param_initial_beta, which_model, all_ratings, all_seq_vals, all_output,nbins_reward,experiment,log_or_not);

ll_beta = 0;

num_subs = size(all_ratings,2);

max_subs = num_subs;

for num_subs_found = 1:max_subs;
    
    %peel off this sub
    seq_vals = all_seq_vals(:,:,num_subs_found);
    mean_ratings = all_ratings(:,num_subs_found);   %returns to original naming convention, at risk of being confusing
    
    %Need for assigning reward
    if log_or_not == 1;
        binEdges_reward = linspace(min(log(mean_ratings(:))),max(log(mean_ratings(:))),nbins_reward+1);   %organise bins by min and max
    else
        binEdges_reward = linspace(min(mean_ratings(:)),max(mean_ratings(:)),nbins_reward+1);   %organise bins by min and max
    end;
    
    %set this sub's prior
    distOptions = 0;
    if log_or_not == 1;
        dataPrior.mean = mean(log(mean_ratings));   %take the log to normalize
        dataPrior.var = var(log(mean_ratings));
        dataPrior.range = [min(min(log(seq_vals))) max(max(log(seq_vals)))];   %used for BV normilisation
    else;
        dataPrior.mean = mean(mean_ratings);   %take the log to normalize
        dataPrior.var = var(mean_ratings);
        dataPrior.range = [min(min(seq_vals)) max(max(seq_vals))];   %used for BV normilisation
    end;
    dataPrior.kappa = 2;
    dataPrior.nu = 1;
    %     options = optimset('Display','iter','MaxFunEvals', 5000, 'TolFun', 0.001);
    
    %now loop through the sequences, just to set up the data to be passed
    %into fit_inner
    for sequence=1:size(seq_vals,1);
        
        %ranks for this sequence
        dataList = tiedrank(seq_vals(sequence,:)');
        
        %Subject draws and ranks
        sub_choice_trial(num_subs_found,sequence) = all_output(sequence,num_subs_found);  %just transferring draws from new variable name back to old one
        sub_choice_rank(num_subs_found,sequence) = dataList(sub_choice_trial(num_subs_found,sequence)); %ranks
        
        %format data for model fitting
        %initialze list
        if log_or_not == 1;
            list(sequence).allVals = log(seq_vals(sequence,:));
            list(sequence).vals = log(seq_vals(sequence,1:sub_choice_trial(num_subs_found,sequence)));
        else
            list(sequence).allVals = seq_vals(sequence,:);
            list(sequence).vals = seq_vals(sequence,1:sub_choice_trial(num_subs_found,sequence));
        end;
        list(sequence).optimize = 0;
        list(sequence).flip = 1;
        list(sequence).length = size( list(sequence).allVals, 2 );
        
    end;    %loop through sequences (will start new sequence loop inside fit_inner to run model on each sequence, this loop here just compiles the data to be passed to fit_inner
    
    if which_model == 1 %cost to sample fit
        
        params = [-0.03];   %initialise Cs & beta
        
    elseif which_model == 2;
        
        %initialise parameter search using psychometric function values
        if experiment == 1; %av
            load('C:\matlab_files\fiance\power_param_study2_6bins_indsub1.mat','parameters');  %pre-fit and used in publication
        elseif experiment == 2; %matchmaker
            load('C:\matlab_files\fiance\power_param_study9_6bins_indsub1.mat','parameters');  %pre-fit and used in publication
        elseif experiment == 3; %trustworthiness
            load('C:\matlab_files\fiance\power_param_study13_6bins_indsub1.mat','parameters');  %pre-fit and used in publication
        end;
        
        params = [parameters(:,num_subs_found)']; %initialise with this subs' data and beta=5;
        
    elseif which_model == 3 | which_model == 4;    %ideal observer and biased values pre-fit will fit only beta (prefit logistic paraneters fixed inside sectretary progtramme)
        
        params = [0]; %I'll just pass in a dummy parameter which does affect computations any so I can use the same function for testing.
        
    elseif which_model == 5;    %biased prior
        
        params = [0.5]; %initialise difference from mean and beta too
        
    elseif which_model == 6;    %cut-off
        
        params = (1/exp(1))*list(1).length; %initialise cut off to 37% rule
        
    end;
    
    %     disp(sprintf('fitting model %d to subject %d', which_model, num_subs_found));
    if num_subs_found == 1;
        fprintf('subject: 1, ');
    elseif num_subs_found == max_subs;
        fprintf('%d ... done\n',num_subs_found);
    else
        fprintf('%d, ',num_subs_found);
    end;
    
    %     options = optimset('Display','iter');
    [mparams(num_subs_found,:), lla(num_subs_found), exitflag(num_subs_found), search_out(num_subs_found)] = ...
        fminsearch(@(params) fit_inner(params, dataPrior, list, distOptions,binEdges_reward,which_model,experiment,num_subs_found,param_initial_beta), params);
    
end;    %loop through subjects

ll_beta = sum(lla);




function analyze_value_position_functions(value_data,choice_trial,plot_cmap,binEdges_psi,legend_labels,param_to_fit,two_params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Look at proportion choice, position and value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

param_to_fit = [0 param_to_fit];    %the first model in this function is subject so add it to this list as zeroth model so it can be indexed (e.g., in colormaps)

%value_data was log(raw_seq_subs), contains ratings data in sequences and is seq*position*sub
nbins = size(binEdges_psi,2)-1;
num_subs = size(value_data,3);
num_positions = size(value_data,2);
num_seqs = size(value_data,1);
num_models = size(choice_trial,3);
serial_r = value_data;  %only for ratings plots, nothing else

f_a = 0.1; %face alpha
sw = 0.5;  %ppoint spread width
font_size = 12;

%ok. this time, let's be less efficient but more organised. I want to bin
%things right up front before anything so there is a ratings dataset
%(value_data) and a binned dataset (value_bins)
for sub = 1:num_subs;
    binEdges = binEdges_psi(sub,:);
    [dummy,value_bins(:,:,sub)] = histc(value_data(:,:,sub), [binEdges(1:end-1) Inf]);
end;

%nan mask has zeros for view decisions and 1s for take decisions and NaNs
%when no face was seen because of elapsed decision. Can use to mask other
%arrays later ...
nan_mask = NaN(num_seqs,num_positions,num_subs,num_models);
for model=1:num_models;
    for sub=1:num_subs;
        for seq=1:num_seqs;
            nan_mask(seq,1:choice_trial(sub,seq,model),sub,model) = 0;
            nan_mask(seq,choice_trial(sub,seq,model),sub,model) = 1;
        end;
    end;
end;

%now we have two seq*position*subject arrays, one ratings, one bins, now
%make supersubject seq*position versions by concatenating subjects. Effectively, each
%new subject justs adds new sequences so its a long list of sequences
value_data_sub = [];
value_bins_sub = [];
choice_trial_sub = [];

for sub=1:num_subs;
    value_data_sub = [value_data_sub; value_data(:,:,sub)];
    value_bins_sub = [value_bins_sub; value_bins(:,:,sub)]; %bins are still subject specific
    choice_trial_sub = [choice_trial_sub; squeeze(choice_trial(sub,:,:))];
end;

nan_mask_sub = [];
for model=1:num_models;
    temp = [];
    for sub=1:num_subs;
        temp = [temp; squeeze(nan_mask(:,:,sub,model))];
    end;
    nan_mask_sub(:,:,model) = temp;
end;

%yes, it's yet another subject loop. I'm being modular. This one prepares
%the average ratings * serial position data. It also computes the proportion choices *serial position data.
%It also computes proportion subject predicted * serial position data
%This one needed whether using a super subject or fitting all subjects, it's a separate analysis
model_predicted_choices = NaN(num_subs,num_positions,num_models-1); %for proportion correctly predicted subject choices
position_choices = NaN(num_subs,num_positions,num_models);   %for proportion responses serial positions
position_function = zeros(num_positions,num_subs,num_models);   %for average ratings as function of serial position
position_it = zeros(num_positions,num_subs,num_models);         %for average ratings as function of serial position
for sub=1:num_subs;
    
    this_subject_ratings = squeeze(serial_r(:,:,sub));  %only for plotting ratings by serial position
    
    for position=1:num_positions; %loop through the positions
        for model=1:num_models;
            
            sub_choices_this_position = nan_mask(:,position,sub,1);
            model_choices_this_position = nan_mask(:,position,sub,model);
            % %             %computes proportion responses for each position
            position_choices(sub,position,model) = sum( choice_trial(sub,:,model) == position )/size(choice_trial,2);
            
            %find average attractiveness of the choices in each position
            this_subject_choices = squeeze(choice_trial(sub,:,model));
            indices_into_choices = find(this_subject_choices==position);
            if ~isempty(indices_into_choices);
                for i=1:size(indices_into_choices,2);
                    position_function(position,sub,model) = position_function(position,sub,model)+this_subject_ratings(indices_into_choices(i),position);
                    position_it(position,sub,model) = position_it(position,sub,model)+1;
                end;    %loop through values for this position
            end;    %is there a valid position value?
            
        end;    %model
    end;        %position
end;            %subject

%individual subs
position_data_indiv = NaN(nbins,num_positions,num_subs,num_models);            %for value function analyses as function of serial position
ave_rating_per_bin_indiv = NaN(nbins,num_positions,num_subs,num_models);       %use this for curve fitting later

for sub=1:num_subs;
    
    for position=1:num_positions; %loop through the positions
        for model=1:num_models;
            
            this_subject_bins = value_bins(:,:,sub);
            temp2 = squeeze(nan_mask(:,:,sub,model));
            this_subject_bins(isnan(temp2(:)))=NaN;
            
            for val_bin = 1:nbins;
                
                %find bins at this position and,if any, check what are the CHOICES and RATINGS for that bin/position
                trials_with_bins_in_this_position = [];
                trials_with_bins_in_this_position = find( this_subject_bins(:,position) == val_bin );   %on which sequences did a value in this bin occur in this position?
                num_trial_with_vals = numel(trials_with_bins_in_this_position);                     %how many sequences have this value in this position?
                position_data_indiv(val_bin,position,sub,model) = sum(choice_trial(sub,trials_with_bins_in_this_position,model)==position)/ num_trial_with_vals ; %Now I need the number of CHOICES for this positon and bin
                ave_rating_per_bin_indiv(val_bin,position,sub,model) = nanmean(value_data(trials_with_bins_in_this_position,position,sub));   %need this for regression later (no sense of model here)
                
            end;    %value bin
        end;    %model
    end;        %position
end;            %subject


position_data = position_data_indiv;
ave_rating_per_bin = ave_rating_per_bin_indiv;

%loop again, this time through positions and fit averages over subjects
for model=1:num_models;
    for position=1:num_positions;
        
        %computes value slopes for each position, and model
        this_position_no_subs = nanmean(squeeze( position_data(:,position,:,model) ),2);   %returns bin values for a subject in a position
        
        y_this_position = this_position_no_subs(~isnan(this_position_no_subs));
        x_this_position = [1:nbins];
        x_this_position = x_this_position(~isnan(this_position_no_subs))';
        clear f_position;
        if numel(x_this_position)<3 | position==num_positions | sum(this_position_no_subs) == 0;    %if there are too many nans and not enough datapoints, if its that last position with the flat line or all ones, or if no response was ever made
            
            b1(position,model) = NaN;
            b2(position,model) = NaN;
        else
            %             f_position=fit(x_this_position,y_this_position,'1./(1+exp(-p1*(x-p2)))','StartPoint',[1 5],'Lower',[0 1],'Upper',[Inf 8]);
            %             temp_coef = coeffvalues(f_position);
            %             b1(position,model) = temp_coef(1);  %if only slope and mid are free
            %             b2(position,model) = temp_coef(2); %if only slope and mid are free
            if two_params == 1;
                %             Two params free
                f_position=fit(x_this_position,y_this_position,'1./(1+exp(-p1*(x-p2)))','StartPoint',[1 5],'Lower',[0 1],'Upper',[Inf 8]);
                temp_coef = coeffvalues(f_position);
                b1(position,model) = temp_coef(1);  %if only slope and mid are free
                b2(position,model) = temp_coef(2); %if only slope and mid are free
                
            else
                
                % %             %Three params free
                f_position=fit(x_this_position,y_this_position,'p1./(1+exp(-p3*(x-p4)))','StartPoint',[1 1 5],'Lower',[0 0 1],'Upper',[1 Inf 8]);
                temp_coef = coeffvalues(f_position);
                b1(position,model) = temp_coef(2);  %if only slope and mid are free
                b2(position,model) = temp_coef(3); %if only slope and mid are free
            end;
            
        end;    %check is there enough data to do a fit
        
    end;    %loop through positions
end;    %loop through models

b_ci = zeros(size(b2));
b = b2;

%%%%%%%new part: correlation position data for each model with subjects
r_graph = zeros(1,num_models);
r_ci_graph = zeros(1,num_models);
for model = 1:num_models;
    
    if model ~=1;   %subject is model 1
        
        for sub=1:num_subs;
            
            clear this_subject_data this_model_data this_subject_data_rs this_model_data_rs
            
            %extract_data
            this_subject_data = squeeze(position_data(:,:,sub,1));
            this_model_data = squeeze(position_data(:,:,sub,model));
            
            %reshape data
            this_subject_data_rs = reshape(this_subject_data,prod(size(this_subject_data)),1);
            this_model_data_rs = reshape(this_model_data,prod(size(this_model_data)),1);
            
            %correlate them
            try
                [temp1 temp2] = corrcoef(this_subject_data,this_model_data,'rows','complete');  
            catch
                disp('insufficient data for subject-model correlation');
                temp1 = NaN(4);
                temp2 = NaN(4);
            end;
            r(sub,model-1) = temp1(2,1);
            p(sub,model-1) = temp2(2,1);
            sub_nums(sub,model-1) = sub;
            mod_nums(sub,model-1) = model-1;
            
        end;    %loop through subs
        
    end;    %only consider models other than subjects
    
end;    %models

r_graph = [0 nanmean(r)];
r_ci_graph = [0 1.96*(std(r)/sqrt(size(r,1)))];


%average proportion responses, over subjects
mean_position_choices = squeeze(mean(position_choices));
ci_position_choices = squeeze(1.96*(std(position_choices)/sqrt(size(position_choices,1))));
%average ratings as function of serial position
clear ave_ratings ave ci;
ave_ratings = position_function./position_it;

%serial postion plots: average rating, proportion correct and value sensitivity slopes
h3 = figure; set(gcf,'Color',[1 1 1]);  %For serial position/PSE plots/correlation plots
h4 = figure; set(gcf,'Color',[1 1 1]);  %For psychometric function plots
for model = 1:size(choice_trial,3);
    
    markersize = 3;
    
    %average rating as function of serial positon
    legend_locs = [0.5:-0.05:(0.5 - (0.05*5))];
    
    %proportion choices
    figure(h3); subplot( 2,2,3); hold on;
    sph = shadedErrorBar(1:size(mean_position_choices,1),mean_position_choices(:,model),ci_position_choices(:,model),{'MarkerFaceColor',plot_cmap(param_to_fit(model)+1,:),'MarkerEdgeColor',plot_cmap(param_to_fit(model)+1,:),'Marker','o','MarkerSize',markersize,'LineStyle','-'},1); hold on;
    set(sph.mainLine,'Color',plot_cmap(param_to_fit(model)+1,:));
    set(sph.patch,'FaceColor',plot_cmap(param_to_fit(model)+1,:));
    set(sph.edge(1),'Color',plot_cmap(param_to_fit(model)+1,:));
    set(sph.edge(2),'Color',plot_cmap(param_to_fit(model)+1,:));
    %         text(3,legend_locs(model),legend_names{model},'Color',plot_cmap(param_to_fit(model)+1,:),'FontSize',12,'FontName','Arial');
    box off;
    %axis square;
    set(gca,'FontSize',12,'FontName','Arial','xtick',[1:size(b,1)],'ytick',[0.1:0.1:0.8],'Ylim',[0 0.5],'Xlim',[1 size(b,1)],'LineWidth',2);
    xlabel('Position in Sequence'); ylabel('Proportion Choices');
    
    %psychometric function parameters
    figure(h3); subplot( 2,1,1 ); hold on;
    errorbar(1:size(b,1),b(:,model),b_ci(:,model),'Color',plot_cmap(param_to_fit(model)+1,:),'MarkerFaceColor',plot_cmap(param_to_fit(model)+1,:),'MarkerEdgeColor',plot_cmap(param_to_fit(model)+1,:),'Marker','o','MarkerSize',markersize,'LineStyle','-','LineWidth',1); hold on;
    box off;
    set(gca,'FontSize',12,'FontName','Arial','xtick',[1:size(b,1)],'Xlim',[1 size(b,1)],'LineWidth',2);
    xlabel('Position in Sequence'); ylabel('Point of Subjective Equality');
    
    if model ~=1;   %no subjects
        %model correlations with proportion choice
        figure(h3); subplot( 2,2,4 ); hold on;
        legend_positions = [1.1:-.05:0];
        
        handles = plotSpread(r(:,model-1), ...
            'xValues',model,'distributionColors',plot_cmap(param_to_fit(model)+1,:),'distributionMarkers','.', 'spreadWidth', sw);
        
        bar(model,r_graph(model), ...
            'FaceColor',plot_cmap(param_to_fit(model)+1,:),'FaceAlpha',f_a,'EdgeColor',[0 0 0] );
        
        text(1.5,legend_positions(model),legend_labels{param_to_fit(model)}, 'Color',plot_cmap(param_to_fit(model)+1,:), 'FontSize',12,'FontName','Arial');
        
        set(gca,'FontSize',12,'FontName','Arial', 'xticklabel',{[]},'LineWidth',2);
        ylabel('Model-participant Correlation');
        ylim([0 1.0]);
        xlim([1 numel(r_graph)+0.5]);
        
    end;    %If not a subject
    
    %value psychometric functions (in a different figure with different colormap), with lines for each position
    figure(h4); subplot(1,num_models,model);
    pm_line_colors = cool(size(position_data,2)+1);
    
    for position_line = 1:size(position_data,2)-1;
        
        if numel(size(position_data))==4;
            h = plot( nanmean(squeeze(position_data(:,position_line,:,model)),2) ); hold on;
        else
            h = plot( position_data(:,position_line,model) ); hold on;
        end;
        axis square;
        set(h,'Marker','o','MarkerSize',6,'MarkerEdgeColor',pm_line_colors(position_line,:),'MarkerFaceColor',pm_line_colors(position_line,:),'Color',pm_line_colors(position_line,:),'LineStyle','-','LineWidth',2);
        set(gca,'FontSize',12,'FontName','Arial','xtick',[1:size(position_data,1)],'xlim',[0.5 size(position_data,1)+0.5],'ylim',[0 1.1],'ytick',[0:0.2:1],'LineWidth',2);
        xlabel('Attractiveness Bin'); ylabel('Proportion Choices'); box off;
        
    end;    %position lines
    
    if model == num_models;
        legend('Position 1','Position 2','Position 3','Position 4','Position 5','Position 6','Position 7','Position 8','Position 9','Position 10','Position 11','Position 12');
    end;
    
end;    %loop through models






%%
function [mean_ratings seq_vals output] = get_sub_data(experiment);

%In this version, we'll extract all the subject data here before we embark
%on the subject loop in the main body. That means another subject loop
%inside of here, unfortunately.

if experiment == 1; %online holidays
    
    %sub num then data (rating) for holidays then foods then female faces than male faces
    temp_r = xlsread('C:\matlab_files\fiance\online_domains_01\online_replication\all_attractiveness_ratings_noheader_online.xlsx');
    %There are some nans hanging on the end of domains with less than the max num of subjects. Get rid of them nice are early (here)
    data = temp_r(~isnan(temp_r(:,1)),1:2);
    %sub num then data (draws) for holidays then foods then female faces than male faces
    [temp_d dummy dummy] = xlsread('C:\matlab_files\fiance\online_domains_01\online_replication\all_sequences_online.xlsx');
    draws = temp_d(~isnan(temp_d(:,1)),1:2);
    %These cols are 1:domain code (1=holidays, 2=food, 3=female, 4=male), 2:event num, 3: sub_id, 4-11: the 8 filenames for this sequence
    [dummy dummy temp_f] = xlsread('C:\matlab_files\fiance\online_domains_01\online_replication\all_domains_online_sequence_fnames.xlsx');
    seq_fnames = temp_f(find(cell2mat(temp_f(:,1)) == 1),:);
    %filename index key (cols: holiday, food, male, female) 90 filenames key to ratings
    [dummy1 dummy2 key_ratings] = xlsread('C:\matlab_files\fiance\online_domains_01\domains_key_ratings_fnames.xlsx');
    key_this_domain = key_ratings(:,1);
    
elseif experiment == 2; %online foods
    
     %sub num then data (rating) for holidays then foods then female faces than male faces
    temp_r = xlsread('C:\matlab_files\fiance\online_domains_01\online_replication\all_attractiveness_ratings_noheader_online.xlsx');
    %There are some nans hanging on the end of domains with less than the max num of subjects. Get rid of them nice are early (here)
    data = temp_r(~isnan(temp_r(:,3)),3:4);
    %sub num then data (draws) for holidays then foods then female faces than male faces
    [temp_d dummy dummy] = xlsread('C:\matlab_files\fiance\online_domains_01\online_replication\all_sequences_online.xlsx');
    draws = temp_d(~isnan(temp_d(:,3)),3:4);
    %These cols are 1:domain code (1=holidays, 2=food, 3=female, 4=male), 2:event num, 3: sub_id, 4-11: the 8 filenames for this sequence
    [dummy dummy temp_f] = xlsread('C:\matlab_files\fiance\online_domains_01\online_replication\all_domains_online_sequence_fnames.xlsx');
    seq_fnames = temp_f(find(cell2mat(temp_f(:,1)) == 2),:);
    %filename index key (cols: holiday, food, male, female) 90 filenames key to ratings
    [dummy1 dummy2 key_ratings] = xlsread('C:\matlab_files\fiance\online_domains_01\domains_key_ratings_fnames.xlsx');
    key_this_domain = key_ratings(:,2);
    
    elseif experiment == 4; %open day holidays
    
    %sub num then data (rating) for holidays then foods then female faces than male faces
    temp_r = xlsread('C:\matlab_files\fiance\online_domains_01\open_day_replication\all_attractiveness_ratings_noheader_openday.xlsx');
    %There are some nans hanging on the end of domains with less than the max num of subjects. Get rid of them nice are early (here)
    data = temp_r(~isnan(temp_r(:,1)),1:2);
    %sub num then data (draws) for holidays then foods then female faces than male faces
    [temp_d dummy dummy] = xlsread('C:\matlab_files\fiance\online_domains_01\open_day_replication\all_sequences_openday.xlsx');
    draws = temp_d(~isnan(temp_d(:,1)),1:2);
    %These cols are 1:domain code (1=holidays, 2=food, 3=female, 4=male), 2:event num, 3: sub_id, 4-11: the 8 filenames for this sequence
    [dummy dummy temp_f] = xlsread('C:\matlab_files\fiance\online_domains_01\open_day_replication\all_domains_openday_sequence_fnames.xlsx');
    seq_fnames = temp_f(find(cell2mat(temp_f(:,1)) == 1),:);
    %filename index key (cols: holiday, food, male, female) 90 filenames key to ratings
    [dummy1 dummy2 key_ratings] = xlsread('C:\matlab_files\fiance\online_domains_01\domains_key_ratings_fnames.xlsx');
    key_this_domain = key_ratings(:,1);
    
elseif experiment == 5; %open day foods
    
    %sub num then data (rating) for holidays then foods then female faces than male faces
    temp_r = xlsread('C:\matlab_files\fiance\online_domains_01\open_day_replication\all_attractiveness_ratings_noheader_openday.xlsx');
    %There are some nans hanging on the end of domains with less than the max num of subjects. Get rid of them nice are early (here)
    data = temp_r(~isnan(temp_r(:,3)),3:4);
    %sub num then data (draws) for holidays then foods then female faces than male faces
    [temp_d dummy dummy] = xlsread('C:\matlab_files\fiance\online_domains_01\open_day_replication\all_sequences_openday.xlsx');
    draws = temp_d(~isnan(temp_d(:,3)),3:4);
    %These cols are 1:domain code (1=holidays, 2=food, 3=female, 4=male), 2:event num, 3: sub_id, 4-11: the 8 filenames for this sequence
    [dummy dummy temp_f] = xlsread('C:\matlab_files\fiance\online_domains_01\open_day_replication\all_domains_openday_sequence_fnames.xlsx');
    seq_fnames = temp_f(find(cell2mat(temp_f(:,1)) == 2),:);
    %filename index key (cols: holiday, food, male, female) 90 filenames key to ratings
    [dummy1 dummy2 key_ratings] = xlsread('C:\matlab_files\fiance\online_domains_01\domains_key_ratings_fnames.xlsx');
    key_this_domain = key_ratings(:,2); 
    
elseif experiment == 3; %faces online, a pain (slightly) because we want to concatenate subjects.
    
    %sub num then data (rating) for holidays then foods then female faces than male faces
    temp_r = xlsread('C:\matlab_files\fiance\online_domains_01\online_replication\all_attractiveness_ratings_noheader_online.xlsx');
    %There are some nans hanging on the end of domains with less than the max num of subjects. Get rid of them nice are early (here)
    data = [temp_r(~isnan(temp_r(:,5)),5:6); temp_r(~isnan(temp_r(:,7)),7:8)];  %this time, concatenate columsn for females, then males
    %sub num then data (draws) for holidays then foods then female faces than male faces
    [temp_d dummy dummy] = xlsread('C:\matlab_files\fiance\online_domains_01\online_replication\all_sequences_online.xlsx');
    draws = [temp_d(~isnan(temp_d(:,5)),5:6); temp_d(~isnan(temp_d(:,7)),7:8)];
    %These cols are 1:domain code (1=holidays, 2=food, 3=female, 4=male), 2:event num, 3: sub_id, 4-11: the 8 filenames for this sequence
    [dummy dummy temp_f] = xlsread('C:\matlab_files\fiance\online_domains_01\online_replication\all_domains_online_sequence_fnames.xlsx');
    seq_fnames = temp_f(find(cell2mat(temp_f(:,1)) == 3 | cell2mat(temp_f(:,1)) == 4),:);
    %list of female faces
    female_faces = zeros(size(seq_fnames,1),1);
    female_faces( 1:numel(find(cell2mat(temp_f(:,1)) == 3) )) = 1;
    %filename index key (cols: holiday, food, female, male) 90 filenames key to ratings
    [dummy1 dummy2 key_ratings] = xlsread('C:\matlab_files\fiance\online_domains_01\domains_key_ratings_fnames.xlsx');
    key_this_domain_f = key_ratings(:,3);   %separate list for male faces 
    key_this_domain_m = key_ratings(:,4);   %separate list for female faces
    
    
elseif experiment == 6; %faces openday, a pain (slightly) because we want to concatenate subjects.
    
    %sub num then data (rating) for holidays then foods then female faces than male faces
    temp_r = xlsread('C:\matlab_files\fiance\online_domains_01\open_day_replication\all_attractiveness_ratings_noheader_openday.xlsx');
    %There are some nans hanging on the end of domains with less than the max num of subjects. Get rid of them nice are early (here)
    data = [temp_r(~isnan(temp_r(:,5)),5:6); temp_r(~isnan(temp_r(:,7)),7:8)];  %this time, concatenate columsn for females, then males
    %sub num then data (draws) for holidays then foods then female faces than male faces
    [temp_d dummy dummy] = xlsread('C:\matlab_files\fiance\online_domains_01\open_day_replication\all_sequences_openday.xlsx');
    draws = [temp_d(~isnan(temp_d(:,5)),5:6); temp_d(~isnan(temp_d(:,7)),7:8)];
    %These cols are 1:domain code (1=holidays, 2=food, 3=female, 4=male), 2:event num, 3: sub_id, 4-11: the 8 filenames for this sequence
    [dummy dummy temp_f] = xlsread('C:\matlab_files\fiance\online_domains_01\open_day_replication\all_domains_openday_sequence_fnames.xlsx');
    seq_fnames = temp_f(find(cell2mat(temp_f(:,1)) == 3 | cell2mat(temp_f(:,1)) == 4),:);
    %list of female faces
    female_faces = zeros(size(seq_fnames,1),1);
    female_faces( 1:numel(find(cell2mat(temp_f(:,1)) == 3) )) = 1;
    %filename index key (cols: holiday, food, female, male) 90 filenames key to ratings
    [dummy1 dummy2 key_ratings] = xlsread('C:\matlab_files\fiance\online_domains_01\domains_key_ratings_fnames.xlsx');
    key_this_domain_f = key_ratings(:,3);   %separate list for male faces 
    key_this_domain_m = key_ratings(:,4);   %separate list for female faces
    
end;    %test which experiment and get correct data


num_seq_pos = 8;
these_sub_ids = unique( cell2mat(seq_fnames(:,3)) );

for sub_id = 1:numel(these_sub_ids);
    
    %We're going to replace seq filenames with ratings. Ready the sequence filenames
    this_sub_indices = find( cell2mat(seq_fnames(:,3)) == these_sub_ids(sub_id) );  %index into one of the subs
    this_sub_data = seq_fnames(this_sub_indices,4:11);  %filenames for sequences for this sub
    this_sub_num_seqs = numel(this_sub_indices);    %Should always be 6 sequences
    if exist('key_this_domain_f');  %if rating filenames were divided into male and female earlier
        sexes_this_subject = female_faces(this_sub_indices,1);
        if sexes_this_subject(1) == 1;
            key_this_domain = key_this_domain_f;
        else
            key_this_domain = key_this_domain_m;
        end;
        
    end;    %check if i need to assign sex specific face key
    
    %We're going to replace seq filenames with ratings. Ready the ratings
    this_sub_ratings_indices = find(data(:,1)==these_sub_ids(sub_id)); %This should find all the ratings for this subject in the same domain as the sequences. Should be 90
    this_sub_ratings = data(this_sub_ratings_indices,2);   %ratings themselves for this subject, should be 90. Should be sorted by key
    
    mean_ratings(:,sub_id) = this_sub_ratings;  %first output of this function
    
    %Get number of samples for this subject's sequence (Will use in sequence loop to get ranks for each sequence)
    clear this_sub_draws;
    this_sub_draws = draws(find(draws(:,1)==these_sub_ids(sub_id)),2);    %Should return the 5/6 sequence indices; subids in draws file better match up with subids in sequience lists (they should)
    
    %next function output
    output(:,sub_id) = this_sub_draws;
    
    for seq_num = 1:this_sub_num_seqs;
        for seq_pos = 1:num_seq_pos;
            
            filename_to_find = seq_fnames(this_sub_indices(seq_num),seq_pos+3); %What is this filename?
            idx = find(strcmp(key_this_domain, filename_to_find));  %find the index corresponding to the sequence position's filename
            this_seq_pos_rating = this_sub_ratings(idx);    %get rating corresponding to the index for this filename
            
            %next main output
            seq_vals(seq_num,seq_pos,sub_id) = this_seq_pos_rating;
            
        end;    %loop through sequence positions
    end;    %loop through sequences
end;   %loop through subs


%%
function ll = fit_inner(params, dataPrior, input_list, distOptions,binEdges_reward,param_to_fit, experiment, num_subs_found,param_initial_beta)

%The function formerly known as f_fitparams

%Now beta is no longer passed in with parameters but as a fized variable
%param_initial_beta. Params will be a zero for ideal observer and empirical
%bias values because they don't estimate any parameters here. It should
%just converge right away but entered parameters don't affect output.

b = param_initial_beta;

nTrials = size(input_list,2);

ll = 0;
% b = 0.5;

for triali = 1 : nTrials;
    
    list.length = input_list(triali).length;
    list.allVals = input_list(triali).allVals;
    list.vals = input_list(triali).vals;
    list.flip = input_list(triali).flip;
    list.optimize = input_list(triali).optimize;
    
    listDraws = numel(list.vals);
    
    if param_to_fit == 6;   %If cut-off model
        
        %put ones for every candidate choice
        this_seq_vals = list.allVals;
        choiceStop = zeros(size(this_seq_vals));
        choiceStop(1,(find( this_seq_vals(ceil(params(1))+1:end) > max(this_seq_vals(1:ceil(params(1))))))+ceil(params(1))) = 1;
        choiceStop(1,list.length) = 1;
        %put zeros for every non-candidate
        choiceCont = ones(size(this_seq_vals));
        choiceCont(1,find(choiceStop==1))=0;
        difval = choiceCont-choiceStop;
        
    else;
        %     [choiceStop, choiceCont] = analyzeSecretary_brunoversiontest_2020(dataPrior, list, params, distOptions,param_to_fit);
        [choiceStop, choiceCont difval] = analyzeSecretaryNick_2020(dataPrior, list, params, distOptions,binEdges_reward,param_to_fit,experiment,num_subs_found);
    end;
    
    choiceValues = [choiceCont' choiceStop'];
    
    cprob = zeros(listDraws, 2);
    
    for drawi = 1 : listDraws
        cprob(drawi, :) = exp(b*choiceValues(drawi, :))./sum(exp(b*choiceValues(drawi, :)));
    end
    
    %ll = ll - sum(log(cprob(1:(listDraws-1), 1))) - log(cprob(listDraws, 2));
    if listDraws == 1;
        ll = ll - 0 - log(cprob(listDraws, 2));
    else
        ll = ll - sum(log(cprob((listDraws-1), 1))) - log(cprob(listDraws, 2));
    end;
    
end;


