function fh = plotScatterX1(condfitX, condsimX)

% make scatter plots of simulated and and fitted parameter values (with line
% of fit and correlation coefficient) for cost to sample for each beta

% -----------

% define initial parameters
ylm     = [0 80];
rows    = 2;
cols    = 5;
c       = [[171,5,32]/256;[12,35,75]/256; [92, 135, 39]/256; [132, 210, 226]/256; [7, 104, 115]/256; [241, 158, 31]/256;
    [183, 85, 39]/256; [74, 48, 39]/256; [25, 88, 120]/256; [102, 12, 48]/256];
clf


%% Plot simulated-estimated Cost to sample 

% loop over beta groups 
groups  = length(condfitX);
cc      = 0;

for i = 1:groups

    cc                  = cc + 1;

    Cs_sim              = condsimX{1,i}(1,:,:); Cs_sim = Cs_sim(:);
    Cs_fit              = condfitX{1,i}(1,:,:); Cs_fit = Cs_fit(:);
    ax(cc)              = subplot(rows,cols,cc);

    % add scatter 
    fh                  = scatter(ax(cc),Cs_sim,Cs_fit,[],c(i,:),'*');

    % add line of fit 
    lfit(cc)            = lsline(ax(cc));
    lfit(cc).Color      = c(i,:);
    lfit(cc).LineWidth  = 0.7;


    % add x,y labels 
    ylabel('estimated Cost-to-Sample value')
    xlabel('simulated Cost-to-Sample value')

end 

end % end of function