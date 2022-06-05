% re-write the event values 

% re-writing event values will help epoching draws as draw choices and urn
% choices for the easy and difficult conditions 

% old event values
% 1. easy blue
% 2. easy green
% 3. diff blue
% 4. diff green

% new event values 
% 1. easy/draw
% 2. easy/urn
% 3. diff/draw
% 4. diff/urn

% extract fields that need to be re-written
load('fdfMspmeeg_sub_01_beads_block_01.mat');

% remove first 2 rows of events struct
D.trials.events(:,1:2) = [];

% extract events from struct
allevents   = D.trials.events;

% difne a few vars
trialend    = 103; % trigger code for end of trial
blocknum    = 1;

% extract event values 
eventvalues = [allevents(:).value]';

% length of events array
eventlength = length(eventvalues);

% re-write the event values
for i = 1:eventlength
    
    % if the event trigger is easy (1 or 2)
    if eventvalues(i) == 1 | eventvalues(i) == 2
        
        if eventvalues(i + 4) == trialend % if it is the last trigger/draw
            
            eventvalues(i) = 2; % urn-choice event trigger
        else
            eventvalues(i) = 1; % draw-choice event trigger
        end
        
    % if the event trigger is difficult (3 or 4)
    elseif eventvalues(i) == 3 | eventvalues(i) == 4
        
        if eventvalues(i + 4) == trialend % if it is the last trigger/draw
            
            eventvalues(i) = 4; % urn-choice event trigger
        else
            eventvalues(i) = 3; % draw-choice event trigger
        end
    end
end

% move the new event values array back to the events struct and add the new
% block field
for j = 1:eventlength
    
    D.trials.events(:,j).value = eventvalues(j,1);
    D.trials.events(:,j).block = blocknum;

end

clear allevents blocknum eventlength eventvalues i j trialend

save fdfMspmeeg_sub_01_beads_block_01.mat
