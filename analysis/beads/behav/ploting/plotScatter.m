function f  = plotScatter(y, x)

% make scatter plots of simulated and and fitted sampling rates (with line
% of fit and correlation coefficient)

% --------------

% define initial parameters
colours = [1 0 1]; % let's start with magenta?
rows    = 2;
cols    = 5;
ylm     = [0 10];
xlm     = ylm;
clf

% loop over groups (that is over beta values 
groups = size(x,1);

for ii = 1:groups

    % plot the easy and diff conditions seperately
    ax(ii)          = subplot(rows,cols,ii);

    yg              = y(ii,:); % get mean across reps/subjects and transpose so that the results is a 1xn array
    xg              = x(ii,:);

    % add scatter 
    f               = scatter(ax(ii),xg,yg,'magenta','*');
    
    % add line of fit 
    lfit(ii)        = lsline(ax(ii));
    lfit(ii).Color  = 'magenta';
    lfit(ii).LineWidth  = 0.7;

    % axis labels
    ylabel('estimated number of samples')
    xlabel('simulated number of samples')
        
end % end of groups (beta param values) loop

end % end of function 