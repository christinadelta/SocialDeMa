function f  = plotScatter2(y, x)

% x = this_sim
% y = this_fit

% -------------

% define initial params
rows    = 10;
cols    = 9;
ylm     = [0 10];

c = [[171,5,32]/256;[12,35,75]/256; [92, 135, 39]/256; [132, 210, 226]/256; [7, 104, 115]/256; [241, 158, 31]/256;
        [183, 85, 39]/256; [74, 48, 39]/256; [25, 88, 120]/256; [102, 12, 48]/256];
% clf

% loop over groups (that is over beta values 
groups = size(x,2);
cc      = 0; % initiate counter

for i = 1:groups

    beta_sim    = x{1,i};
    beta_fit    = y{1,i};

    % loop over cs groups (that is over cs values 
    gcs         = size(beta_sim,1);

    for ii = 1:gcs

        cc              = cc + 1; % update counter

        % plot the easy and diff conditions seperately
        ax(cc)          = subplot(rows,cols,cc);

        yg              = beta_fit(ii,:); % get samples for a given comb of beta + cs for all reps/subs
        xg              = beta_sim(ii,:);

        % add scatter plot
        f               = scatter(ax(cc),xg,yg, [],c(i,:),'*');
        

        % add line of fit 
        lfit(cc)            = lsline(ax(cc));
        lfit(cc).Color      = c(i,:);
        lfit(cc).LineWidth  = 0.7;


        % axis labels
%         ylabel('estimated samples')
%         xlabel('simulated samples')

    end % end of cs loop
end % end of betas loop

end % end of function 