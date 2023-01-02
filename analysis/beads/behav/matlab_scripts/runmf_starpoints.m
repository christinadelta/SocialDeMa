function [minParams, ll, Qsad, cprob] = runmf_starpoints(R, cond, thiscond_seqmat, thiscond_choiceVec)

% this function model fitting using two free parameters [cost-to-sample and
% softmax beta for choice probabilities] using different beta values to
% compare fit 

% fit free parameter using Bruno's version of the model
[minParams, ll, Qsad, cprob]        = bayesbeads_b(thiscond_seqmat, thiscond_choiceVec, R);


end