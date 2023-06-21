function fGrid = plotGrid(conditions, outmat, allCs, allbetas)

% plot 2d landscape of of NLL for all combinations of Beta and Cs

fGrid{1}    = figure('Name','Likelihoods easy condition');
set(fGrid{1},'position',[10 60 650 650],'paperunits','centimeters','Color','w');

fGrid{2}    = figure('Name','Likelihoods difficult condition');
set(fGrid{2},'position',[10 60 650 650],'paperunits','centimeters','Color','w');

% loop over conditions
for i = 1:conditions

    % extract condition mat 
    tmpOut = outmat{1,i};

    figure(fGrid{i});
    dims    = conditions;
    %subplot(dims, dims, 1)
    imagesc(tmpOut); hold on;

    % find the best ll
    % [maxCs maxBeta] = find(tmpOut==min(tmpOut(:)));
    [maxCs maxBeta] = find(tmpOut==min(tmpOut(tmpOut>1)));
    hline(maxCs, {'w-', 'LineWidth', 2});
    vline(maxBeta, {'w-', 'LineWidth', 2});

    title('p(data|model)')
    ytick       = allCs;
    xtick       = allbetas;
    xticklabel  = num2cell(xtick(2:end),1);
    yticklabel  = num2cell(ytick,1);
    set(gca,'ytick',ytick,'yticklabel',yticklabel,'xtick',xtick,'xticklabel',xticklabel)
    ylabel('Cost sample'); xlabel('beta')


end % end of conditions loop

end 

