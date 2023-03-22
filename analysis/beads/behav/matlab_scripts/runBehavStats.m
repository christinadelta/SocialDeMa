function [draws_stats,acc_stats, pc_results, pc_tables] = runBehavStats(nsubs,all_draws,all_acc,all_ioacc,all_iodraws)

% created: 21/03/2023

% INPUTS: 
%           - number of subjects (1x1 double)
%           - all_draws (2xnsbus) array with averaged draws for each subject and conditon
%           - all_acc (2xnsbus) array with averaged performance for each subject and conditon
%           - all_iodraws (2xnsbus) array with averaged ideal-observer draws 
%           - all_ioacc (2xnsbus) array with averaged ideal-observer performance

% THE FUNCTION RUNS:
% 1) RM ANOVAs using probabilitiy and agent type as factors
% 2) pairwise comparisons 


% OUTPUT:
%           - RM ANOVA output
%           - Pairwise comparisons output

% -------

% define variables

nconds              = size(all_acc,2);
totaltrials         = 52;
condtrials          = 26;


% make arrays to be used in the ANOVAS
subvec              = repmat(1:nsubs,1,4)';                             % create a vector with 4 copies participant number 
agentvec            = repmat([ones(1,nsubs*2) ones(1,nsubs*2)*2],1,1)'; % create a vector with 2 copies of agent type (indexed as 1=human, 2=io)
probvec             = repmat([ones(1,nsubs) ones(1,nsubs)*2],1,2)';     % create a vector with 2 copies of probability type (indexed as 1=0.8, 2=0.6)

% create 1 vec with all draws (human, io) 
drawsmat            = [all_draws all_iodraws];
drawsvec            = drawsmat(:);

% create 1 vec with all acc (human, io) 
accmat              = [all_acc all_ioacc];
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


end % end of function