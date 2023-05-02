function [output_struct] = runBehavStats(nsubs, anova_struct)

% created: 21/03/2023

% INPUTS: 
%           - number of subjects (1x1 double)
%           - structure with agent data:     
%           - all_draws (2xnsbus) array with averaged draws for each agent type and conditon
%           - all_acc (2xnsbus) array with averaged performance for each agent type and conditon

% THE FUNCTION RUNS:
% 1) RM ANOVAs using probabilitiy and agent type as factors
% 2) pairwise comparisons 


% OUTPUT:
%           - RM ANOVA output
%           - Pairwise comparisons output

% -------

% what is the size of the input structure? needed to know the levels of agent type factor
struct_size = length(fieldnames(anova_struct));

% define variables
% nconds              = size(all_humanacc,2);
totaltrials         = 52;
condtrials          = 26;


if struct_size <= 4
    

    % unpack structure 
    all_humandraws = anova_struct.all_draws; all_iodraws = anova_struct.all_iodraws; 
    all_humanacc = anova_struct.all_acc; all_ioacc = anova_struct.all_ioacc;

    % make arrays to be used in the ANOVAS
    subvec              = repmat(1:nsubs,1,4)';                             % create a vector with 4 copies participant number 
    agentvec            = repmat([ones(1,nsubs*2) ones(1,nsubs*2)*2],1,1)'; % create a vector with 2 copies of agent type (indexed as 1=human, 2=io)
    probvec             = repmat([ones(1,nsubs) ones(1,nsubs)*2],1,2)';     % create a vector with 2 copies of probability type (indexed as 1=0.8, 2=0.6)

    % create 1 vec with all draws (human, io) 
    drawsmat            = [all_humandraws all_iodraws];
    drawsvec            = drawsmat(:);
    

    % create 1 vec with all acc (human, io) 
    accmat              = [all_humanacc all_ioacc];
    accvec              = accmat(:);

    % make tables vor vis 
    acctable            = table(subvec,agentvec,probvec,accvec);
    drawstable          = table(subvec,agentvec,probvec,drawsvec);

    %% RUN BEHAV STATISTICS - RM ANOVAs %%

    % run rm 2x2 anova on draws 
    [pvals,~,stats]     = anovan(drawsvec, {subvec agentvec probvec}, ... 
    'model','interaction', 'random',1,'varnames',{'subvec' 'agentvec' 'probvec'});
    
    % store anova results 
    draws_stats.pvals   = pvals;
    draws_stats.stats   = stats;
    
    clear pvals stats
    
    % run rm 2x2 anova on accuracy 
    [pvals,~,stats]     = anovan(accvec, {subvec agentvec probvec}, ... 
    'model','interaction', 'random',1,'varnames',{'subvec' 'agentvec' 'probvec'});
    
    % store anova results 
    acc_stats.pvals     = pvals;
    acc_stats.stats     = stats;

    %% PAIRWISE COMPARISONS %%

    % 1. multicompare draws
    % run multicompare on agent_type and probability factors
    [draws_results,~,~,groups]  = multcompare(draws_stats.stats,"Dimension",[2 3]); 
    % [c, m, h, nms] = multcompare(draws_stats.stats,"Dimension",[2 3],'Alpha',.05, 'CriticalValueType','bonferroni'); 
    
    % look at results in a table 
    draws_tbl                   = array2table(draws_results,"VariableNames",...
        ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
    
    % change groups to group-names 
    draws_tbl.("Group A")       = groups(draws_tbl.("Group A"));
    draws_tbl.("Group B")       = groups(draws_tbl.("Group B"));
    
    % 2. multicompare accuracy 
    % run multicompare on agent_type and probability factors
    [acc_results,~,~,groups]    = multcompare(acc_stats.stats,"Dimension",[2 3]); 
    
    % look at results in a table 
    acc_tbl                     = array2table(acc_results,"VariableNames",...
        ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
    
    % change groups to group-names 
    acc_tbl.("Group A")         = groups(acc_tbl.("Group A"));
    acc_tbl.("Group B")         = groups(acc_tbl.("Group B"));
    
    % store output
    pc_results.draws            = draws_results;
    pc_results.acc              = acc_results;
    pc_tables.draws             = draws_tbl;
    pc_tables.acc               = acc_tbl;

    output_struct = struct('draws_stats',draws_stats,'acc_stats',acc_stats,...
        'pc_results',pc_results,'pc_tables',pc_tables);

elseif struct_size >= 4

    % unpack structure 
    all_humandraws      = anova_struct.all_draws; 
    all_iodraws         = anova_struct.all_iodraws; 
    all_model1_draws    = anova_struct.costSample_modelSamples;
    all_model2_draws    = anova_struct.beta_modelSamples;
    all_model3_draws    = anova_struct.betaCs_modelSamples;

    all_humanacc        = anova_struct.all_acc; 
    all_ioacc           = anova_struct.all_ioacc;
    all_model1_acc      = anova_struct.costSample_modelPerf;
    all_model2_acc      = anova_struct.beta_modelPerf;
    all_model3_acc      = anova_struct.betaCs_modelPerf;

    % make arrays to be used in the ANOVAS
    subvec              = repmat(1:nsubs,1,10)';                             % create a vector with 10 copies participant number 
    agentvec            = repmat([ones(1,nsubs*2) ones(1,nsubs*2)*2 ones(1,nsubs*2)*3 ones(1,nsubs*2)*4 ...
        ones(1,nsubs*2)*5],1,1)'; % create a vector with 2 copies of agent type (indexed as 1=human, 2=io, 3=beta, 4=beta_Cs)
    probvec             = repmat([ones(1,nsubs) ones(1,nsubs)*2],1,5)';     % create a vector with 2 copies of probability type (indexed as 1=0.8, 2=0.6)
    
    % create 1 vec with all draws (human, io and parameterised models) 
    drawsmat            = [all_humandraws all_iodraws all_model1_draws all_model2_draws all_model3_draws all_model4_draws all_model5_draws all_model6_draws];
    drawsvec            = drawsmat(:);
    
    % create 1 vec with all acc (human, io and parameterised models) 
    accmat              = [all_humanacc all_ioacc all_model1_acc all_model2_acc all_model3_acc all_model4_acc all_model5_acc all_model6_acc];
    accvec              = accmat(:);

    % make tables vor vis 
    acctable            = table(subvec,agentvec,probvec,accvec);
    drawstable          = table(subvec,agentvec,probvec,drawsvec);

    %% RUN BEHAV STATISTICS - RM ANOVAs %%

    % run rm 2x2 anova on draws 
    [pvals,~,stats]     = anovan(drawsvec, {subvec agentvec probvec}, ... 
    'model','interaction', 'random',1,'varnames',{'subvec' 'agentvec' 'probvec'});

    % store anova results 
    draws_stats.pvals   = pvals;
    draws_stats.stats   = stats;
    
    clear pvals stats

    % run rm 2x2 anova on accuracy 
    [pvals,~,stats]     = anovan(accvec, {subvec agentvec probvec}, ... 
    'model','interaction', 'random',1,'varnames',{'subvec' 'agentvec' 'probvec'});
    
    % store anova results 
    acc_stats.pvals     = pvals;
    acc_stats.stats     = stats;

    %% PAIRWISE COMPARISONS %%

    % 1. multicompare draws
    % run multicompare on agent_type and probability factors
    [draws_results,~,~,groups]  = multcompare(draws_stats.stats,"Dimension",[2 3]); 
    
    
    % look at results in a table 
    draws_tbl                   = array2table(draws_results,"VariableNames",...
        ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
    
    % change groups to group-names 
    draws_tbl.("Group A")       = groups(draws_tbl.("Group A"));
    draws_tbl.("Group B")       = groups(draws_tbl.("Group B"));

    % 2. multicompare accuracy 
    % run multicompare on agent_type and probability factors
    [acc_results,~,~,groups]    = multcompare(acc_stats.stats,"Dimension",[2 3]); 
    
    % look at results in a table 
    acc_tbl                     = array2table(acc_results,"VariableNames",...
        ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
    
    % change groups to group-names 
    acc_tbl.("Group A")         = groups(acc_tbl.("Group A"));
    acc_tbl.("Group B")         = groups(acc_tbl.("Group B"));

    % store output
    pc_results.draws            = draws_results;
    pc_results.acc              = acc_results;
    pc_tables.draws             = draws_tbl;
    pc_tables.acc               = acc_tbl;

    output_struct = struct('draws_stats',draws_stats,'acc_stats',acc_stats,...
        'pc_results',pc_results,'pc_tables',pc_tables);

end % end of if statement 


end % end of function