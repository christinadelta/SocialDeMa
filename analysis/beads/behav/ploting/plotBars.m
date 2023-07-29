function h = plotBars(mdl_fitsamples)

% created 10/07/2023
% plots bar graphs with individual points


%% plot for each condition

rows=1;
cols=2;
ylm = [0 10];
clf

for jj = 1:2 % 2 conditions 

    % plot the easy condition first
    subplot(rows,cols,jj)

    mdl_easy    = mdl_fitsamples(:,:,jj);
    mdl_easyAV  = mean(mdl_easy,2)'; % get mean across reps/subjects and transpose so that the results is a 1xn array
    h(jj)        = bar(mdl_easyAV, 0.4, 'hist');
    set(h(jj),'FaceColor',[1 0 1]) % change colour of bars to green
    
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
    % xticklabels('beta1', 'beta2','beta3','beta4','beta5','beta6','beta7','beta8','beta9','beta10')

end % end of conditions loop


end % end of function 