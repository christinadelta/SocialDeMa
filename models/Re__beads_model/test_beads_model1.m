%%% TEST SCRIPT FOR NICK'S BEADS MODEL 

% create simulated data 
sequence1   = [1 1 0 1 0 0 1 1 0 1];
sequence2   = [1 1 1 0 0 1 1 0 1 0];

seq_mat     = [sequence1; sequence2];
% seq_mat = [sequence1];
% seq_mat = [sequence2];


% need sequence and behavior, but behavior needs to be in its old format
alpha           = 1;        % softmax stochasticity parameter (for fitting to human behaviour)
Cw              = -1000;    % The difference between the rewards for being correct (in this case no reward 0) and the cost of being wrong (-1000).
q               = 0.6;      % proportion of the majority value in sequence (60/40 split in this case)
Cs              = -10;      % the cost to sample
aqvec_switch    = 1;        % still not sure why exactly this is needed 