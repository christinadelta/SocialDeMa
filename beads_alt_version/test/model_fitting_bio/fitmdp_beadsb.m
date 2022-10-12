function [ll] = fitmdp_beadsb(param, R, seq_mat, choiceVec, condraws)

[ll, Qsad, cprob] = fitmdp_beads(param, R, seq_mat, choiceVec, condraws);


return