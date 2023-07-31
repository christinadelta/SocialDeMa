function gh = plotScatterX2(condfitX, condsimX)

% make scatter plots of simulated and and fitted parameter values (with line
% of fit and correlation coefficient)

%%% NOT USED RIGHT NOW 

% -----------

% define initial parameters
ylm     = [0 80];
rows    = 1;
cols    = 2;
c       = [[171,5,32]/256;[12,35,75]/256; [92, 135, 39]/256; [132, 210, 226]/256; [7, 104, 115]/256; [241, 158, 31]/256;
    [183, 85, 39]/256; [74, 48, 39]/256; [25, 88, 120]/256; [102, 12, 48]/256];
clf


%% Plot simulated-estimated Betas 

cc      = 0;
% loop over beta groups 
groups  = length(condfitX);

for i = 1:groups

    cc                      = cc + 1;

    beta_sim                = condsimX{1,i}(2,:,:); beta_sim = beta_sim(:);
    beta_fit                = condfitX{1,i}(2,:,:); beta_fit = beta_fit(:);
    beta_sim_all(:,i)       = beta_sim; 
    beta_fit_all(:,i)       = beta_fit; 
   
end

% flatten x and y
beta_sim_all = beta_sim_all(:);
beta_fit_all = beta_fit_all(:);

ax                  = subplot(rows,cols,1);


% add scatter 
gh                      = scatter(ax,beta_sim_all,beta_fit_all,[],c(1,:),'*');

% add line of fit 
lfit            = lsline(ax);
lfit.Color      = c(1,:);
lfit.LineWidth  = 0.7;

% rescale the y axis 
if ~any(isnan(ylm(1,:))), ylim(ylm(1,:)); end
ym      = get(gca,'ylim');


% add x,y labels 
ylabel('estimated beta value')
xlabel('simulated beta value')

 
end % end of function