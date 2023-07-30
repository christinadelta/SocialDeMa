function fh = plotScatterX(y, x)

% make scatter plots of simulated and and fitted parameter values (with line
% of fit and correlation coefficient)

% --------

% define initial parameters
colours = [1 0 1]; % let's start with magenta?
rows    = 1;
cols    = 2;
clf

% loop over conditions
for jj = 1:2

    mdl_fitx        = y(:,:,jj); mdl_fitx = mdl_fitx(:); % flatten
    mdl_simx        = x(:,:,jj); mdl_simx = mdl_simx(:);

    ax(jj)          = subplot(rows,cols,jj);

    % add scatter 
    fh              = scatter(ax(jj),mdl_simx,mdl_fitx,'magenta','*');

    % add line of fit 
    lfit(jj)            = lsline(ax(jj));
    lfit(jj).Color      = 'magenta';
    lfit(jj).LineWidth  = 0.7;

    % add x,y labels 
    ylabel('estimated beta value')
    xlabel('simulated beta value')

end

end % end of function