function  makePlots(beta, cs, ll, draws)

subplot(2,2,1)
plot(beta)
title('optimal beta values')
xlabel('subjects')
ylabel('estimated beta')
legend('easy (0.8)', 'difficult (0.6)')

subplot(2,2,2)
plot(cs)
title('optimal cost-to-sample values')
xlabel('subjects')
ylabel('estimated cost-sample')
legend('easy (0.8)', 'difficult (0.6)')

subplot(2,2,3)
plot(ll)
title('minimum negative ll values')
xlabel('subjects')
ylabel('estimated negative ll')
legend('easy (0.8)', 'difficult (0.6)')

subplot(2,2,4)
plot(draws)
title('subject average draws')
xlabel('subjects')
ylabel('draws')
legend('easy (0.8)', 'difficult (0.6)')


return