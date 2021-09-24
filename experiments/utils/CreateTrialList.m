function [trials, set] = CreateTrialList(set)

% The function uses the task number to create a list of trials for each task 
% part of the OPTIMAL STOPPING experiments

% extract task number 
taskNb                      = set.taskNb;

if taskNb == 1
    
    % UNPACK THE SETTINGS STRUCTURE
    total_draws                 = set.draws; 
    diff_conds                  = set.conds;
    total_trials                = set.trials;
    trials_cond                 = total_trials / diff_conds;
    blocks                      = set.blocks;
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

    % create sequences of draws for the two conditions 
    for cond = 1:diff_conds
        for i = 1:trials_cond

            this_prob               = ceil(p(cond) * total_draws); % if p =0.8, then high_prob = 8
            tempseq                 = cat(2, ones(1,this_prob), ones(1,total_draws - this_prob)*2);

            if length(triallist) >= trials_cond

                i                   = length(triallist) + 1; % increment i so that it doesn't overwrite the trials.sequence struct
                triallist{i}        = tempseq(randperm(total_draws));

            else
                triallist{i}        = tempseq(randperm(total_draws));

            end % end of if statement
        end % end of sequence loop
    end % end of conditions loop

    % shuffle the sequences 
    order                           = randperm(numel(triallist));
    triallist                       = triallist(order);

    temp                            = 0; % for splitting trials in blocks 
    % split sequences and urns in blocks 
    for block = 1:blocks

        trials.sequence{block}      = triallist(1 + temp:blocktrials*block);
        trials.urns{block}          = urns(1 + temp:blocktrials*block);
        temp                        = temp + blocktrials; % update

    end % end of blocks loop
    
elseif taskNb == 2 % if this is the economic task
    
    % which phase is it?
    phase = set.phase;
    
    if phase == 1 % if this is the ratings phase 
        
        % UNPACK SETTINGS
        totaltrials     = set.totaltrials;
        blocks          = set.blocks;
        blocktrials     = set.blocktrials; 
        itemReps        = set.itemReps;
        
        tempArray      = 1:totaltrials/2; % create an array with all the trials
        trialArray     = repmat(tempArray,1,itemReps);
        
        order           = randperm(numel(trialArray));% randomise the trial array
        trialArray      = trialArray(order); % 
        
        temp            = 0; % for splitting trials in blocks 
        
        % trials split in blocks 
        for i = 1:blocks
            
            trials.sequence{i}      = trialArray(1 + temp:blocktrials*i);
            temp                    = temp + blocktrials; % update i
            
            
        end
        
    else % if phase == 2
        
        % UNPACK SETTINGS
        blocks          = set.blocks;       % number of blocks
        samples         = set.samples;      % number of samples per trial
        totaltrials     = set.totaltrials;       % number of trials
        blocktrials     = set.blocktrials;  % number of trials per block
        phaseitems      = set.phaseitems;   % number of items/contracts in the 2nd phase (65% of phase 1 items)
        items           = set.items;        % total contracts
        
        % choose a 65% subset from the total contracts/items (that should
        % be 300)
        arraysize       = numel(items);
        idx             = randperm(arraysize);
        templist        = items(idx(1:phaseitems));
        
        % split the list in 30 sequences of 10 samples 
        temp            = 0; 
        for s = 1:totaltrials 
            
            list{s}     = templist(1 + temp:samples*s);
            temp        = temp + samples;
            
        end
        
        % no split sequences in blocks (10 per block)
        temp            = 0; % for splitting trials in blocks 
        for i = 1:blocks
            
            trials.sequence{i}      = list(1 + temp:blocktrials*i);
            temp                    = temp + blocktrials;
            
        end
   
    end % end of phase if statement 
    
elseif taskNb == 3
    
    % which phase is it?
    phase = set.phase;
    
    if phase == 1
        
        % UNPACK SETTINGS
        totaltrials     = set.totaltrials;
        blocks          = set.blocks;
        blocktrials     = set.blocktrials; 
        itemReps        = set.itemReps;
        
        tempArray       = set.items'; % create an array with all the trials
        trialArray      = repmat(tempArray,1,itemReps);
        
        order           = randperm(numel(trialArray));% randomise the trial array
        trialArray      = trialArray(order); % 
        
        temp            = 0; % for splitting trials in blocks 
        
        % trials split in blocks 
        for i = 1:blocks
            
            trials.sequence{i}      = trialArray(1 + temp:blocktrials*i);
            temp                    = temp + blocktrials; % update i
            
            
        end
        
    else % if this is phase 2
        
        % UNPACK SETTINGS
        blocks          = set.blocks;       % number of blocks
        samples         = set.samples;      % number of samples per trial
        totaltrials     = set.totaltrials;  % number of trials
        blocktrials     = set.blocktrials;  % number of trials per block
        phaseitems      = set.phaseitems;   % number of items in the 2nd phase (62% of phase 1 items)
        items           = set.items(:,1);        % total faces
        
        % choose a 62% subset from the total contracts/items (that should
        % be 300)
        arraysize       = numel(items);
        idx             = randperm(arraysize);
        templist        = items(idx(1:phaseitems));
        
        % split the list in 30 sequences of 10 samples 
        temp            = 0; 
        for s = 1:totaltrials 
            
            list{s}     = templist(1 + temp:samples*s);
            temp        = temp + samples;
            
        end
        
        % no split sequences in blocks (10 per block)
        temp            = 0; % for splitting trials in blocks 
        for i = 1:blocks
            
            trials.sequence{i}      = list(1 + temp:blocktrials*i);
            temp                    = temp + blocktrials;
            
        end
        
    end % end of phase statement 
    
end % end of taskNb if statement 

end