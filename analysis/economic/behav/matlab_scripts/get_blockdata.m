function [sequences,choicevecs,phase2_data, samples] = get_blockdata(subdir,sub,task)

% This function is used with Phase 1 -- Economic data

% the function runs through economic_prepro_v3.m and onwards
% created in December 2022
% extracts data from logs and stores in all_data matrix and subsequences
% cell for further processing & analysis

session     = 2;
phase       = 2;
counter     = 0; % init counter var
blocks      = 2;
blocktrials = 20;
respoptions = 2; % accept vs decline (for phase 2)

for block = 1:blocks

    fprintf('\t\t loading block %d\n\n',block);
    subFile = fullfile(subdir, sprintf('subject_%02d_task_%s_block_%02d_ses_%02d_phase_%02d_blocktrials_logs.mat',sub,task,block,session,phase));
    load(subFile)

    for trial = 1:blocktrials
            
        indx                                = ((block -1)*blocktrials) + trial; 
        id                                  = ((block -1)*blocktrials) + trial;
        
        subj(indx)                          = sub;
        blockno(indx)                       = block;
        trialno(indx)                       = logs.blocktrials(trial).trialnumber;
        numsamples(indx)                    = logs.blocktrials(trial).numsamples;
        thisitem(indx)                      = logs.blocktrials(trial).chosenitem;
        thisprice(indx)                     = logs.blocktrials(trial).chosenprice;
        thisrank(indx)                      = logs.blocktrials(trial).rank;
        reward(indx)                        = logs.blocktrials(trial).reward;
        balance(indx)                       = logs.blocktrials(trial).balance;
        
        sequences{1,id}                     = logs.blocktrials(trial).sequence;
        
        % create a vector of sequence responses at this point. This
        % will be based on the number of samples on every
        % trial/sequence
        t                                   = nan(numsamples(indx), respoptions); % init empty vec
        
        for s = 1:numsamples(indx) 
            
            if s < numsamples(indx) % if this is a decline response
                t(s,1)                      = 1;
                t(s,2)                      = 0;
            else                    % if participant accepted an option
                t(s,1)                      = 0;
                t(s,2)                      = 1;
            end
                       
        end % end of samples loop
        
        % store temporal vec (t) in cell
        choicevecs{1,id}                    = t;
        
        clear t s 
        
    end % end of trials loop


end % end of blocks loop

phase2_data = [subj' blockno' trialno' numsamples' thisitem' thisprice' thisrank' reward' balance'];

% % save matrix in csv format for r and python
% csvwrite('economic_phase2_blockdata.csv', phase2_blockdata)

clear subj trialno blockno thisitem thisprice numsamples thisrank indx

samples = phase2_data(:,4);

end