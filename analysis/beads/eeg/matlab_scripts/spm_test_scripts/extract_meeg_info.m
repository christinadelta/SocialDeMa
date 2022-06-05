% useful commands to extract data from the D object 

% to load the object into workspace
Dobject = spm_eeg_load; 

% look at channel labels
chanlabels = D.chanlabels(:);

% look at the list of conditions 
conds = condlist(D);

% extract all events from D object
block_events = conditions(D);

% what is the length of the block events?
tmp = length(block_events);


