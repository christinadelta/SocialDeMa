function [subsequences,subchoiceVec,all_data] = get_blockdata(subdir,subI,task)

% this function runs through beads_prepro_v5.m
% created in October 2022
% extracts blocktrial data from logs and stores in all_data matrix and subsequences
% cell for further processing & analysis

sv              = 1; % this will be used for spliting sequences in conditions one and two 
cv              = 1; % this will be used for spliting sequencesin conditions one and two 
respoptions     = 3; % b,g,s - needed when recovering choiceVecs
maxdraws        = 10;
session         = 1;
blocktrials     = 13;
blocks          = 4;

for blockI = 1:blocks
    
    fprintf('\t\t loading block %d\n\n',blockI);
    subFile = fullfile(subdir, sprintf('subject_%02d_task_%s_block_%02d_ses_%02d_logs.mat',subI, task, blockI,session));
    load(subFile)
    
    % exctrat trial data
    for trial = 1:blocktrials
        
        indx                            = ((blockI -1)*blocktrials) + trial; 

        block(indx)                     = blockI;
        trialno(indx)                   = logs.blocktrials(trial).trialnumber;
        urntype(indx)                   = logs.blocktrials(trial).urntype;
        
        % here, unfortunately we don't include the last draw (which is
        % the urn choice), the draw_count variable in the experiment
        % was storing only when subject was pressing draw-again...
        if logs.blocktrials(trial).draws == maxdraws
            draws(indx)                 = logs.blocktrials(trial).draws - 1; % this is because of a mistake when coding the experiment (if it was the last bead and sub pressed draw-again, this variable was updated! shit! 
        elseif logs.blocktrials(trial).draws < maxdraws
            draws(indx)                 = logs.blocktrials(trial).draws;
        end
        response(indx)                  = logs.blocktrials(trial).response;
        accuracy(indx)                  = logs.blocktrials(trial).accuracy;
        condition(indx)                 = logs.blocktrials(trial).condition;
        subj(indx)                      = subI;
        generaltrl(indx)                = indx;

        sequence                        = logs.blocktrials(trial).sequence; % extract tish trial sequence of draws
        
        % maybe here add this trial's responses. Given that I was not saving within sequence resposes but we need that info for
        % the model-fitting part, I will create a nx3 matrix of choices for each sequence, when n=number of draws+choice and 3=choice options (b,g,d)
        draws_choice                    = draws(indx)+1; % + 1 for choice (very last draw)

        t                               = nan(draws_choice, respoptions); % init empty matrix
        
        for d = 1:draws_choice
            if d ~= draws_choice % if this is not the last draw add 0's to b and g columns and 1 to draw column
                t(d,1:2)                = 0; % index zero for b and g columns
                t(d,3)                  = 1; % index one for draw column
            else
                if urntype(indx) == 1 & accuracy(indx) == 1 % this is a blue urn and sub was correct
                    t(d,2:3)            = 0; % index zero for g and draw columns
                    t(d,1)              = 1; % index one for b column (sub ressed blue)
                elseif urntype(indx) == 1 & accuracy(indx) == 0 % this is a blue urn and sub was incorrect or did not respond
                    t(d,1)              = 0; % index zero for b 
                    t(d,2)              = 1; % index one for g column
                    t(d,3)              = 0; % index zero for draw 
                elseif urntype(indx) == 0 & accuracy(indx) == 1 % this is a green urn and sub was correct
                    t(d,1)              = 0; % index zero for b 
                    t(d,2)              = 1; % index one for g column
                    t(d,3)              = 0; % index zero for s 
                elseif urntype(indx) == 0 & accuracy(indx) == 0 % this is a green urn and sub was incorrect or did not respond
                    t(d,2:3)            = 0; % index zero for g and draw columns
                    t(d,1)              = 1; % index one for b column            
                end
            end % end of if statement
        end % end of draws loop
        
        % add thistrial sequence and choive vector t in correct cell
        % based on condition
        
        if condition(indx) == 1 
            
            subsequences{1,condition(indx)}{1,sv}   = sequence;
            subchoiceVec{1,condition(indx)}{1,sv}   = t;

            sv                                      = sv+1; % update sv
        else
            subsequences{1,condition(indx)}{1,cv}   = sequence;
            subchoiceVec{1,condition(indx)}{1,cv}   = t;

            cv                                      = cv+1; % update cv
  
        end % end of if statement

        clear t sequence  
        
    end % end of trials loop
    
end % end of block trials 

% store all blocks data in a matrix and return
all_data = [subj' block' trialno' urntype' draws' response' accuracy' condition' generaltrl'];

    
return