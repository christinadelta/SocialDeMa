function [trials, set] = CreateTrialList(set)

% The function uses the task number to create a list of trials for each task 
% part of the OPTIMAL STOPPING experiments

% extract task number 
taskNb                      = set.taskNb;

% UNPACK THE SETTINGS STRUCTURE
total_draws                 = set.draws; 
diff_conds                  = set.conds;
total_trials                = set.trials;
trials_cond                 = total_trials / diff_conds;
nb_blocks                   = set.blocks;
blocktrials                 = set.blocktrials;
p                           = set.prob;

set.condtrials              = trials_cond;

% init sequence struct
triallist                   = [];

% create the urns for each trial/sequence 
templist                    = mod(randperm(total_trials),2);
master_order                = randperm(numel(templist)); 
urntemp                     = templist(master_order);
urns                        = urntemp;

%  create sequences of draws for the two conditions 
for cond = 1:diff_conds
    for i = 1:trials_cond
        
        this_prob               = ceil(p(cond) * total_draws); % if p =0.8, then high_prob = 8
        tempseq                 = cat(2, ones(1,this_prob), ones(1,total_draws - this_prob)*2)';
        
        % Add a loss index, when 0 = no loss, when 1 = £10 loss
        if i <=13 
            tempseq(:,2)        = 0;
        else 
            tempseq(:,2)        = 1;
        end 
        
        if length(triallist) >= trials_cond
            
            i                   = length(triallist) + 1; % increment i so that it doesn't overwrite the trials.sequence struct
            triallist{i}        = tempseq(randperm(total_draws),:);
            
        else
            triallist{i}        = tempseq(randperm(total_draws),:);
            
        end % end of if statement
    end % end of sequence loop
end % end of conditions loop

% shuffle the sequences 
order                           = randperm(numel(triallist));
triallist                       = triallist(order);

temp                            = 0; % for splitting trials in blocks 
% split sequences and urns in blocks 
for block = 1:nb_blocks
    
    trials.sequence{block}      = triallist(1 + temp:blocktrials*block);
    trials.urns{block}          = urns(1 + temp:blocktrials*block);
    temp                        = temp + blocktrials; % update
   
end % end of blocks loop

end