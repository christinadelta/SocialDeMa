function hfig = plotCorrSamples (condRfitSamples, condSimSamples, condFitSamples)

% plot samples
hfig = figure;set(hfig,'position',[10 60 500 500],'Color','w'); box off;
[foo,foo,h(1)] = myScatter(condRfitSamples,condSimSamples,false,[0 0 1],'x');
[foo,foo,h(2)] = myScatter(condRfitSamples,condFitSamples,false,[1 0 0],'o');
h(3) = plot([0 1],[0 1],'k:','linewidth',2);
legend(h(1:2),{'Simulated samples','Fit Samples'},'location','best');legend boxoff
% xlim([2 6]); ylim([2 6]);
xlabel('Samples - true'); ylabel ('Samples - estimated');
title('true vs. estimated sampling rates');



end % end of function