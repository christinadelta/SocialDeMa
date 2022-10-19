% Testing fminsearch 

% this is to be used for modelling/fitting the beads data, with teo free
% params

% fminsearch is used to find a minima of a function(s)
% first with 1 variable
% then with two variables


% first find local minima for function f(x)
function f = Ex12_2a_ftn(x)

f = 10 - 5*exp(-(x-0.7).^2) + x.^2;

return

