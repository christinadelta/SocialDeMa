function [allfrontal, allpar] = eegExtract_v3(sub_eeg, sub_drawinfo,sub_AQdiffs, sub_cond)

% VERSION 3 contains code to extract EEG only for parietal and frontal
% responses


% this function is part of the BEADS formal analysis
% Inputs: - a vector that is nx2 (with subjects draws and trial numbers)
%         - structure with cropped EEG data
%         - subject number

% first input (sub_eeg structure) contains data (channels x samples x
% epochs), conditions and events
% we need all these to seperate data for different channels (sites and laterality) and average by epochs/draws

% since the 2nd input is the number of draws, we want to compare with EEG
% data to ensure that behavioural AQ values of the model and EEG data are
% of the same length

% Output: 2 structures based on scalp sites (frontal, parietal) 

% -------------------------------------------
data            = sub_eeg.cropped_struct.data; 
events          = sub_eeg.cropped_struct.events; 
eeg_trials      = length(events);
model_trials    = length(sub_AQdiffs);
condtrials      = 26;


% separete EEG data based on channels 
fdata           = data([1,2,3,4,11,12,13,14,15],:,:);       % all frontal channles
pdata           = data([5,6,7,8,9,10,16,17,18,19,20],:,:);  % all parietal channles

% average the 3d matrices by epochs
av_fdata        = mean(reshape(fdata, [], size(fdata,3)))';
av_pdata        = mean(reshape(pdata, [], size(pdata,3)))';

%% split EEG values to conditions based on events

e               = 1; % easy index
d               = 1; % difficult index 

this_f          = av_fdata;
this_f(:,2)     = events;
this_f(:,3)     = sub_drawinfo(:,1);
this_p          = av_pdata;
this_p(:,2)     = events;
this_p(:,3)     = sub_drawinfo(:,1);

for trl = 1:eeg_trials

    % extract epochs and add them to easy_epochs/dificult cell
    thiseeg_f               = this_f(trl,:);     % frontal site
    thiseeg_p               = this_p(trl,:);     % parietal site

    % split EEG trials into conditions (given that the order of the AQ
    % difference values is: easy trials and then difficult trials)

    % store this-sequence events and eeg in correct condition cell
    if thiseeg_f(1,2) == 1 || thiseeg_f(1,2) == 2 % if this is an easy trial
        fcell{1}(e,:)        = thiseeg_f(:,:); % keep only eeg values, we don't need the events anymore
        pcell{1}(e,:)        = thiseeg_p(:,:);
        
        e                    = e + 1; % update e

    elseif thiseeg_f(1,2) == 3 || thiseeg_f(1,2) == 4 % if this is a difficult trial

        fcell{2}(d,:)        = thiseeg_f(:,:); % keep only eeg values, we don't need the events anymore
        pcell{2}(d,:)        = thiseeg_p(:,:);
        
        d                    = d + 1; % updated d
    end % end of if statements
end % end of trials

%% check if aqdiffs and eeg have same length %%


for cond = 1:2 % conditions


    counter = 0;
    
    cond_front  = fcell{cond,1};
    cond_par    = pcell{cond,1};
    cond_draws  = sub_cond{1,cond}(:,5)+1;
    tmp_AQs     = find(sub_AQdiffs(:,3) == 1);
    cond_AQs    = sub_AQdiffs((tmp_AQs),:);

    new_cond_AQs = nan(size(cond_front,1),1);

    for s = 1:condtrials
        
        % how many draws for this trial?
        this_draws  = cond_draws(s);
        tmpfront    = cond_front(counter+1:counter+this_draws,:);
        tmppar      = cond_par(counter+1:counter+this_draws,:);
        
        % how many AQ values for this trial?
        tmp_diffs   = find(cond_AQs(:,2)==s);
        s_diffs     = cond_AQs((tmp_diffs),:);
        len_diffs   = size(s_diffs,1);
        
        % make up for the difference in AQ values and eeg values
        if len_diffs < size(tmpfront,1)

            comp_diff                                   = size(tmpfront,1) - len_diffs;
            s_diffs(len_diffs+1:len_diffs+comp_diff,:)  = nan;
            ss                                          = size(s_diffs,1);
            new_cond_AQs(counter+1:counter+ss,1)        = s_diffs(:,1);

        elseif len_diffs > size(tmpfront,1)

            % what is the difference?
            newdiff                                     = s_diffs(1:size(tmpfront,1))';
            ss                                          = size(newdiff,1);
            new_cond_AQs(counter+1:ss,1)                = newdiff;

        end % end of if statement

        counter                                         = counter + size(new_cond_AQs,1); % update counter 
        
        clear ss s_diffs newdiff comp_diff 
    end


end 




end % end of function 