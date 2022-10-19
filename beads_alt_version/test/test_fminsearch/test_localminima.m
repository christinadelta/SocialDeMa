% run fmninsearch using this script

clear all 
clc

% Part 1) minimising f = 10 - 5*exp(-(x-0.7).^2) + x.^2
% run part1 of the script to find the minimum f of x 
x = linspace(0,1);
f = Ex12_2a_ftn(x); 
figure 
plot(x,f, 'k', 'Linewidth', 2);
set(gca, 'Fontsize', 24)
set(gcf, 'Paperpositionmode', 'auto')
print(gcf, 'Ex12_2a_graph.jpg', '-djpeg', '-r300')


% CALL FMINSEARCH to find the minimum of x0
% what is the paramater that you want to estimate (find the local minima?)
x0 = 0.35;                   % initial value of the parameter that we want to estimate
fhandle = @Ex12_2a_ftn;     % function to be handled (function that estimates x)
[x, fval] = fminsearch(fhandle,x0); % function handle and initial value are required
disp('minimum of fvalue is')
disp(fval)
disp('... and is located at x = ')
disp(x)
disp('')

% we can use extra params when calling fminsearch
% e.g., x = fmisearch(fhandle, x0, options)
% or x = fmisearch(fhandle, x0, [], param1, param2...)
% output x (where is the minimum) and fval(what is the value of that minimum)

clc
clear all

% Part 2) call fminsearch with more than 1 parameter (param1, param2....)
% what is the value of the minimum x and where is is located?
% CALL FMNISEARCH to find the minimum of Y0
xa  = 10;
ya  = 10;
A   = 10;
    
Yo  = [10;10];
[Y, fval] = fminsearch(@Ex12_2b_ftn,Yo,[],xa,ya, A); 
disp('minimum fvalue is')
disp(fval)
disp('... and is located at x,y = ')
disp(Y)
disp('')

%% lets try part 2 with options included


clc 
clear all

% Set options to monitor the process as fminsearch attempts to locate a
% minimum and to plot the objective function at each iteration.
options = optimset('PlotFcns',@optimplotfval);

xa  = 10;
ya  = 10;
A   = 10;
    
Yo  = [10;10];
[Y, fval] = fminsearch(@Ex12_2b_ftn,Yo,options,xa,ya, A); 
disp('minimum fvalue is')
disp(fval)
disp('... and is located at x,y = ')
disp(Y)
disp('')







