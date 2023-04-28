function [cprob_samples,model_urnchoice] = computeModelSamples(cprob,N)

% compute model sampling rate using cprob instead
% created 24/03/2023

% ---------------------------

len     = 10;                   % true sequence length
trials  = length(cprob(:,1,1)); % number of sequences

% if len 

% loop over sequences 
for t = 1:trials

    choiceProbs             = squeeze(cprob(t,:,:)); % extract this trial choice probabilities 

    % find competing urn - (i.e., the urn with the highest probability to be chosen)
    modelc                  = find(squeeze(cprob(t,:,3)) - max(squeeze(cprob(t,:,1:2))') <0); % which options in this trial have urn > than draw?
    pickdraw                = modelc(1);                            % pick the first one 
    [val urnchoice]         = max(squeeze(cprob(t,pickdraw,:)));    % which of the two urns was chosen?
    model_urnchoice(t,1)    = urnchoice;                            % to be used in the ANOVA analyses

    % sometimes the model draws les than the actual sequence length; find rows with sum zero...
    for l = 1:size(choiceProbs,1)

        if choiceProbs(l,:) == sum(0)
            choiceProbs(l,:) = nan;
        end


    end % end of cprob length loop

    % ...and remove them
    choiceProbs(any(isnan(choiceProbs), 2), :) = [];

    % in cases that the model uses less beads than the total sequence (e.g., 6 draws instaed of 10), the
    % sequence vector we get is 6 rows long; I fill the rest of column 3
    % with -inf and with 1s the column of the competing/chosen urn 
    if size(choiceProbs,1) < len

        if urnchoice == 1 % if model chose the blue urn
            tmp_diff                        = [ones(len,1) zeros(len,1) zeros(len,1)]; % [1 0 0]
        elseif urnchoice == 2 % if model chose the green urn
            tmp_diff                        = [zeros(len,1) ones(len,1) zeros(len,1)]; % [0 1 0]
        end

        tmp_diff(1:size(choiceProbs,1),:)   = choiceProbs;
        choiceProbs                         = tmp_diff;
    end 
    
    % force stop on last draw 
    choiceProbs(end,urnchoice) = Inf;
    
    % ok now its time to compute model's sampling rates based on the 
    for i = 1:N

        isample             = rand(1,size(choiceProbs,1))';

        % now that we have this information let's see what the model is choosing based on the "softmax" choice probabilities 
        tmp_samples(i,1)    = find(choiceProbs(:,urnchoice) > isample(:,1),1, "first");
        
        clear isample
    end % end of sampling loop

    cprob_samples(t,1)      = mean(tmp_samples);

    clear choiceProbs tmp_samples

end % end of trials loop

end % end of function