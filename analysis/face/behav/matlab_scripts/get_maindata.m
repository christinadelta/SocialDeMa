function [subsequences,subchoiceVec,all_data] =  get_maindata(subdir,sub)

% extract phase two block data for the facial attractiveness task
% created 8/2/2023

blocks      = 2;
session     = 2;
phase       = 2;
taskname    = 'face';

% loop over blocks
for blockI = 1:blocks

    fprintf('\t\t loading block %d\n\n',blockI);
    subFile = fullfile(subdir, sprintf('subject_%02d_task_%s_block_%02d_ses_%02d_phase_%02d_blocktrials_logs.mat',sub,taskname,blockI,session,phase));
    load(subFile)

    blocktrials         = length(logs.blocktrials); % how many sequences?

    % loop over trials/sequences
    for trial = 1:blocktrials
            
        indx                                = counter + ((blockI -1)*blocktrials) + trial; 
        id                                  = ((blockI -1)*blocktrials) + trial;
        
        subj(indx)                          = subI;
        blockno(indx)                       = logs.trials(trial).block;
        trialno(indx)                       = logs.trials(trial).trialnumber;
        numsamples(indx)                    = logs.blocktrials(trial).numsamples;
        thisitem(indx)                      = logs.blocktrials(trial).chosenitem;
        
        allsubs_sequences{1,subI}{1,id}     = logs.blocktrials(trial).sequence;
        
    end % end of trials loop


end













end % end of function