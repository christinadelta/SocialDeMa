sequence1 = [1 1 0 1 0 0 1 1 0 1];
sequence2 = [1 1 1 0 0 1 1 0 1 0];
sequence3 = [1 0 0 1 0 1 1 0 1 1];
sequence4 = [0 0 1 1 0 1 1 0 1 1];

seq_mat = [sequence1; sequence2; sequence3; sequence4];
% seq_mat = [sequence1];
% seq_mat = [sequence2];


 %need sequence and behavior, but behavior needs to be in its old format
alpha = 1;  %softmax
Cw = -1000;
q = 0.6;
Cs = -10;

[ll, pickTrial, dQvec, ddec, aQvec, choice] = estimateLikelihoodftest(alpha,Cw,q,Cs,seq_mat,1);

choice(find(choice==2)) = 0;
all_accuracy = mean(choice==1)
all_draws = mean(pickTrial)
all_points = 0 + (sum(choice==0)*-1000) - sum(pickTrial)