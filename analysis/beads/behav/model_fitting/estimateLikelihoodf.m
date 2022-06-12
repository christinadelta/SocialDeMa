function [ll, pickTrial, dQvec, ddec, aQvec] = estimateLikelihoodf(params, sequence, setdata, fixedparams, findpick)

%%% extract parameters from parameter vector
Cw = params(1);
% Cs = params(2);
% alpha  = params(3);
% q = fixedparams;

alpha   = fixedparams(1);
q       = fixedparams(2);
Cs      = fixedparams(3);

lseq    = size(sequence, 2);    % length of sequence
ntrials = size(setdata, 1);     % number of sequences 
ll      = 0;                    % init log likelihood to 0








return