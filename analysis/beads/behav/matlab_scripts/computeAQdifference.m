function AQdiffs = computeAQdifference(sub_model, condtrials, sub)

% this function is part of the BEADS task analyses 
% it computes difference in AQ values for each subject and for each condition
% and checks whether AQ values and participant draws are of the same length

% INPUT:     - moldel output (we need the Q values)
%            - number of condition trials
%            - subject number of draws

% OUTPUT: AQ difference values

% ------------

conds               = 2;

for cond = 1:conds
    
    counter         = 0;
    % extract q values 
    cond_Qs         = sub_model(cond).Qsad; 

    % loop over trials 
    for trl = 1:condtrials

        trialQs     = squeeze(cond_Qs(trl,:,:)); 

        % find competing urn - (i.e., the urn with the highest)
        tmpQs                   = find(squeeze(cond_Qs(trl,:,3)) - max(squeeze(cond_Qs(trl,:,1:2))') <0); % which options in this trial have urn > than draw?
        pickfirst               = tmpQs(1);                            % pick the first one 
        picklast                = tmpQs(end);
        [val urnchoice]         = max(squeeze(cond_Qs(trl,pickfirst,:)));    % which of the two urns was chosen?

        differences                                             = trialQs(:,3) - trialQs(:,urnchoice);

        trialQ_diff_first{cond}(counter+1:counter+pickfirst,1)  = differences(1:pickfirst);
        trialQ_diff_first{cond}(counter+1:counter+pickfirst,2)  = trl;
        trialQ_diff_first{cond}(counter+1:counter+pickfirst,3)  = cond;
        counter                                                 = counter + pickfirst; % update counter 

        
    end % end of trials loop

end % end of conditions loop

AQdiffs = [trialQ_diff_first{1,1}; trialQ_diff_first{1,2}];

% clear trialQ_diff_first

end % end of function