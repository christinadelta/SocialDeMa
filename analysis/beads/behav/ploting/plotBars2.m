function h = plotBars2(this_fit)

% created 26/07/2023
% plots bar graphs of sampling rates with individual points for each beta
% parameter value and Cs value

% -----------

rows = 2;
cols = 5;
ylm = [0 10];

colours = [[171,5,32]/256;[12,35,75]/256; [92, 135, 39]/256; [132, 210, 226]/256; [7, 104, 115]/256; [241, 158, 31]/256;
    [183, 85, 39]/256; [74, 48, 39]/256; [25, 88, 120]/256; [102, 12, 48]/256];
% clf

% how many betas?
groups = length(this_fit);

for jj = 1:groups

    % plot the easy condition first
    subplot(rows,cols,jj)
    
    % extract data for this beta
    mdl_easy    = this_fit{1,jj};
    mdl_easyAV  = mean(mdl_easy,2)'; 

    h(jj)        = bar(mdl_easyAV, 0.4, 'hist');
    set(h(jj),'FaceColor',colours(jj,:)) % change colour of bars to green

    % add individual points 
    hold on
    mdl_ip      = mean(mdl_easy,1); % get individual points 
    
    % plot ips
    for i = 1:length(mdl_easyAV)
        p(i)    = plot(zeros(size(mdl_ip))+i, mdl_easy(i,:),'k.');
    end
    
    % add some jitter to the scatter
    for ii = 1:length(p)
	    x       = get(p(ii),'XData');
	    set(p(ii),'XData',x+randn(size(x))*0.05);
    end
    
    % rescale the y axis 
    if ~any(isnan(ylm(1,:))), ylim(ylm(1,:)); end
    ym      = get(gca,'ylim');

    ylabel('number of samples')


end % end of betas loops

end % end of function 