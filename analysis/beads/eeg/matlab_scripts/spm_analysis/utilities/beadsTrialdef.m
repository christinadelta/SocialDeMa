function [trl, conditionlabels, S] = beadsTrialdef(S)
% Definition of trials based on events for beads task analysed with spm 
% Based on the spm_eeg_trialdef.m function 

% The function is modified to re-write event codes for epoching as:
% - 1. easy draw
% - 2. easy urn
% - 3. difficult draw
% - 4. difficult urn

% format:
% [trl, conditionlabels, S] = beads_definetrial(S)

% INPUT:
% S             - structure

% Fields of S:
% S.D           - filename of EEG mat-file
% S.timewin     - peristimulus time
% S.trialdef    - structure array for trial definition with fields:
%   S.trialdef.conditionlabel - string label for the condition
%   S.trialdef.eventtype      - string (e.g., STATUS)
%   S.trialdef.eventvalue     - number, corressponds to tigger code
%   S.trialdef.trlshift       - shift the triggers by a fixed amount {ms}
% S.save                      - save trial definition in .mat file [1/0]

% OUTPUT:
% trl                        - Nx3 matrix [trlstart trlend offset]
% conditionlabels            - 1xN cell array of strings (labels of events)
% S                          - modified S structure 

% -------------------------------------------------------------------------

% load the meeg object
D = S.D;
obj = spm_eeg_load(D);

% to define trials/epochs, we need continuous data
if ~isequal(obj.type, 'continuous')
    error('trial definition requires continuous dataset as input.');
end

% define parameters 
event       = events(obj, 1, 'samples');
fsample     = obj.fsample;
trialend    = 103;

% define peristimulus time
pretrig     = S.timewin(1);
posttrig    = S.timewin(2);

% calculate trialshift (if trlshift is not 0)
for i = 1:length(S.trialdef)
    if ~isfield(S.trialdef(i),'trlshift')
        trlshift(i) = 0;
    else
        trlshift(i) = round(S.trialdef(i).trlshift * fsample/1000); % assume passed as ms
    end
end

% ectract event values 
% eventvalues = event(:).value';

% calculate the difference in eventvalues and newvalues columns
% diff_values = length(event) - length(eventvalues);

% add (diff_values) columns
% eventvalues = cat(1, zeros(diff_values,1), eventvalues);
newvalues = zeros(length(event),1);


% add new event values to the events struct
for i = 1:length(event)
    
    % if the event trigger is easy (1 or 2)
    if ([event(i).value]) == 1 | ([event(i).value]) == 2
        
        if ([event(i+4).value]) == trialend % if it is the last trigger/draw
            
            newvalues(i) = 2; % urn-choice event trigger
        else
            newvalues(i) = 1; % draw-choice event trigger
        end
        
    % if the event trigger is difficult (3 or 4)
    elseif ([event(i).value]) == 3 | ([event(i).value]) == 4
        
        if ([event(i+4).value]) == trialend % if it is the last trigger/draw
            
            newvalues(i) = 4; % urn-choice event trigger
        else
            newvalues(i) = 3; % draw-choice event trigger
        end
        
    else
        newvalues(i) = 200;

    end
    
    % add new new values as a new field to events struct
    event(:,i).newvalue = newvalues(i);
end


% build trl based on selected events
trl = [];
conditionlabels = {};

for i = 1:numel(S.trialdef)
    
    sel = []; % events selection
    % select all events with specified [i] trigger value and type 
    for j=find(strcmp(S.trialdef(i).eventtype, {event.type}))
        if isempty(S.trialdef(i).eventvalue)
            sel = [sel j];
            
            % Find the values common to both A and B.
        elseif ~isempty(intersect(event(j).newvalue, S.trialdef(i).eventvalue))
            sel = [sel j];
        end
    end
    
    for k=1:length(sel)
        % override the offset of the event
        trloff = round(0.001*pretrig*fsample);        
        % also shift the begin sample with the specified amount
        if ismember(event(sel(k)).type, {'trial', 'average'})
            % In case of trial events treat the 0 time point as time of the
            % event rather than the beginning of the trial 
            trlbeg = event(sel(k)).sample - event(sel(k)).offset + trloff;
        else
            trlbeg = event(sel(k)).sample + trloff;
        end
        trldur = round(0.001*(-pretrig+posttrig)*fsample);
        trlend = trlbeg + trldur;
        
        % Added by Rik in case wish to shift triggers (e.g, due to a delay
        % between trigger and visual/auditory stimulus reaching subject).
        trlbeg = trlbeg + trlshift(i);
        trlend = trlend + trlshift(i);
        
        % add the beginsample, endsample and offset of this trial to the list
        trl = [trl; trlbeg trlend trloff];
        conditionlabels{end+1} = S.trialdef(i).conditionlabel;
    end
      
end % end of trial definition for loop

% sort the trl in right temporal order
[junk, sortind] = sort(trl(:,1));
trl             = trl(sortind, :);
conditionlabels = conditionlabels(sortind);

% modify S structure
S = [];
S.D = D;
S.trl = trl;
S.conditionlabels = conditionlabels;

%%% delete values 200
% newtmp = newvalues 
% newtmp(newtmp==200)=[]

return
