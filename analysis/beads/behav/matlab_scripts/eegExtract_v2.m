function [allfrontal, allfc, allcp, allpar] = eegExtract_v2(sub_eeg, sub_drawinfo)

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

% Output: 4 structures based on scalp sites (frontal, fc, cp, parietal) 
% all 4 structs contain: all (sites) averaged data, left and right (sites)
% averaged data. These will be used to run regressions (with AQ difference as predictor)

% ---------------------------------------------

data            = sub_eeg.cropped_struct.data; 
events          = sub_eeg.cropped_struct.events; 
trials          = length(events);

% separete EEG data based on channels 
% laterality will be: all(1), left(2), right(3). Thus we'll create 4 (sites
% structures with 3 fields

fdata{1}        = data([1,2,3,4,18,19,20,21,22],:,:);           % all frontal channles
fdata{2}        = data([1,2,3,4,18],:,:);                       % left frontal channels
fdata{3}        = data(19:22,:,:);                              % right frontal channels
    
fcdata{1}       = data([5,6,7,23,24,25,26],:,:);                % all frontocentral channels
fcdata{2}       = data([5,6,7,26],:,:);                         % left frontocentral channels
fcdata{3}       = data([23,24,25],:,:);                         % right frontocentral channels

cpdata{1}       = data([8,9,10,17,27,28,29],:,:);               % allcentroparietal channels
cpdata{2}       = data([8,9,10,17],:,:);                        % left centroparietal channels
cpdata{3}       = data([27,28,29],:,:);                         % right centroparietal channels

pdata{1}        = data([11,12,13,14,16,30,31,32,33],:,:); % all parietal channels
pdata{2}        = data([11,12,13,14,16],:,:);                              % left parietal channels
pdata{3}        = data(30:33,:,:);                              % right parietal channels

% pdata{1}        = data([11,12,13,14,15,16,30,31,32,33,34],:,:); % all parietal channels
% pdata{2}        = data(11:16,:,:);                              % left parietal channels
% pdata{3}        = data(30:34,:,:);                              % right parietal channels

%% average the 3d matricies by epochs (3rd dimension)

% for all sites, all, left, right sites average eeg channels and samples for each trial/epoch

% loop over laterality (all, left, right)
for l = 1:3
    
    % first average frontal site
    tmpf            = fdata{l};
    av_fdata{l}     = mean(reshape(tmpf, [], size(tmpf,3)))';
    
    % average fc site
    tmpfc           = fcdata{l};
    av_fcdata{l}    = mean(reshape(tmpfc, [], size(tmpfc,3)))';
    
    % average cp site
    tmpcp           = cpdata{l};
    av_cpdata{l}    = mean(reshape(tmpcp, [], size(tmpcp,3)))';
    
    % average parietal sites
    tmppar          = pdata{l};
    av_pdata{l}     = mean(reshape(tmppar, [], size(tmppar,3)))';
    
    clear tmpf tmpfc tmpcp tmppar
end

%% split EEG values to conditions based on events 

% loop over trials, check if all eeg epochs (for each trial) are of the
% same length as the draws and split into conditions (given that the order of the AQ
% difference values is: easy trials and then difficult trials)

% first loop over laterality (1=all, 2=left, 3=right)
for l = 1:3
    
    % initialise index counters
    e               = 1; % easy index
    d               = 1; % difficult index

    % exctract l array from each of the 4 sites cells and add events
    this_f          = av_fdata{l};
    this_f(:,2)     = events;
    this_fc         = av_fcdata{l};
    this_fc(:,2)    = events;
    this_cp         = av_cpdata{l};
    this_cp(:,2)    = events;
    this_p          = av_pdata{l};
    this_p(:,2)     = events;
    
    for trl = 1:trials
        
        % extract epochs and add them to easy_epochs/dificult cell
        thiseeg_f               = this_f(trl,:);     % frontal site
        thiseeg_fc              = this_fc(trl,:);    % frontocentral site
        thiseeg_cp              = this_cp(trl,:);    % centro-parietal site
        thiseeg_p               = this_p(trl,:);     % parietal site
        
        % store this-sequence events and eeg in correct condition cell
        if thiseeg_f(1,2) == 1 || thiseeg_f(1,2) == 2 % if this is an easy trial
            
            fcell{l}{1}(e,1)        = thiseeg_f(:,1); % keep only eeg values, we don't need the events anymore
            fccell{l}{1}(e,1)       = thiseeg_fc(:,1);
            cpcell{l}{1}(e,1)       = thiseeg_cp(:,1);
            pcell{l}{1}(e,1)        = thiseeg_p(:,1);
            
            e                       = e + 1; % update e
            
        elseif thiseeg_f(1,2) == 3 || thiseeg_f(1,2) == 4 % if this is a difficult trial
            
            fcell{l}{2}(d,1)        = thiseeg_f(:,1); % keep only eeg values, we don't need the events anymore
            fccell{l}{2}(d,1)       = thiseeg_fc(:,1);
            cpcell{l}{2}(d,1)       = thiseeg_cp(:,1);
            pcell{l}{2}(d,1)        = thiseeg_p(:,1);
            
            d                       = d + 1; % updated d
        end % end of if statement
        
    end % end trials loop 
end % end of laterality loop

%% add all (lateralitity, condition, sequences) epochs in the correct order

% now that wehave split the data based on laterlity, conditions and
% sequences for each site and laterality concatinate sequence epochs (first
% easy and then difficult)

% loop over laterality 
for l = 1:3
    
    % extract data from cells
    this_f          = fcell{l};
    this_fc         = fccell{l};
    this_cp         = cpcell{l};
    this_p          = pcell{l};
    
    % exctract condition cells one and two (easy and difficult) and
    % concatinate them
    % frontal
    f_one           = this_f{1};
    f_two           = this_f{2};
    allfrontal{l}   = cat(1,f_one,f_two); % frontal electrodes
    
    % frontocentral
    fc_one          = this_fc{1};
    fc_two          = this_fc{2};
    allfc{l}        = cat(1,fc_one,fc_two); % fc electrodes
    
    % centroparietal
    cp_one          = this_cp{1};
    cp_two          = this_cp{2};
    allcp{l}        = cat(1,cp_one,cp_two); % cp electrodes
    
    % parietal
    p_one           = this_p{1};
    p_two           = this_p{2};
    allpar{l}         = cat(1,p_one,p_two); % p electrodes
    
    clear f_one f_two fc_one fc_two cp_one cp_two p_one p_two
    
end % end of laterality loop

 
end