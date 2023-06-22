function hfig = plotCorrSamples (condRfitSamples, condSimSamples, condFitSamples)

% plot samples
hfig = figure;set(hfig,'position',[10 60 500 500],'Color','w'); box off;
[foo,foo,h(1)] = myScatter(condRfitSamples,condSimSamples,false,[0 0 1],'x');
[foo,foo,h(2)] = myScatter(condRfitSamples,condFitSamples,false,[1 0 0],'o');
h(3) = plot([0 1],[0 1],'k:','linewidth',2);
legend(h(1:2),{'Maximum Likelihood','Expected Value'},'location','best');legend boxoff
xlim([0 1]); ylim([0 1]);
xlabel('p(blue) - true'); ylabel ('p(blue) - estimated');
title('true vs. estimated choice probability');



end % end of function