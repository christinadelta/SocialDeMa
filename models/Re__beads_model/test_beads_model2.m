%% Beads model test Nick's code (Bruno's Version)

sequence1 = [1 1 0 1 0 0 1 1 0 1];
sequence2 = [1 1 1 0 0 1 1 0 1 0];

% define parameters
R.correct = 0;
R.error = -1000;
R.sample = -10;
R.q = 0.6;

seq_mat = [sequence1; sequence2];

% let's define a few variables 
Ntrials         = size(seq_mat,1);
maxDraws        = size(seq_mat,2);
drawSequence    = seq_mat;

%% Run simulated data

[r, Qsa] = backWardInduction(Ntrials, maxDraws, drawSequence, R)





