% test ideal observer (parameter free)

%This code operates on these two sequences, 1s are the positions of the majority outcome
sequence1 = [1 1 0 1 0 0 1 1 0 1];
sequence2 = [1 1 1 0 0 1 1 0 1 0];
sequence3 = [1 0 1 1 0 1 0 0 1 1];
sequence4 = [0 1 1 1 0 0 1 1 0 1];

seq_mat = [sequence1; sequence2; sequence3; sequence4];

%need sequence and behavior, but behavior needs to be in its old format
alpha = 1;   % softmax stochasticity parameter (for fitting to human behaviour)
Cw = -1000;  % The difference between the rewards for being correct (in this case no reward 0) and the cost of being wrong (-1000).
q = 0.6;     % proportion of the majority value in sequence (60:40 split in this case)
Cs = -10;    % the cost to sample

% init arrays to store the models outputs
choices     = zeros(1, size(seq_mat,1));
accs        = zeros(1, size(seq_mat,1));
draws       = zeros(1, size(seq_mat,1));
points      = zeros(1, size(seq_mat,1));


for this_sequence = 1:size(seq_mat,1)
    
    sequence = seq_mat(this_sequence,:);
    
    [ll, pickTrial, dQvec, ddec, aQvec choice] = estimateLikelihoodf(alpha,Cw,q,Cs,sequence,1);
    
    % update model outputs for this sequence
    choices(this_sequence) = choice;
    draws(this_sequence) = pickTrial;

end 

choices(find(choices==2)) = 0;
all_accuracy = mean(choices==1);
all_draws = mean(draws);
all_points = 0 + (sum(choices==0)*Cw) - sum(draws);


% choice(find(choice==2)) = 0;
% all_accuracy = mean(choice==1)
% all_draws = mean(pickTrial)
% all_points = 0 + (sum(choice==0)*-1000) - sum(pickTrial)