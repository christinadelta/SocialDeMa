% run simulated data?

Ntrials  = 4;
maxDraws = 10;
allTrials = 1; % use this to generate simulated sequences 

R.correct = 1;
R.error   = 0;
R.sample  = -0.025;
R.q       = 0.6;
q         = R.q;

for t = 1:Ntrials
    drawSequence = generateDrawSequences(R.q, maxDraws, allTrials);
    allsequences(t,:) = drawSequence;
end

% run backward induction? what is backward induction? calculates utilities
% and action values 
[r, Qsa] = backWardInduction(Ntrials, maxDraws, allsequences, R);

figure;

for dri = 1 : 4
    subplot(2,2,dri);
    plot(1:maxDraws, squeeze(Qsa(dri, :, :)), 'LineWidth', 2);

    choiceTrial = find(squeeze(Qsa(dri, :, 3)) - max(squeeze(Qsa(dri, :, 1:2))') < 0);
    text(choiceTrial(1), 0.1, 'X', 'horizontalalignment', 'center');

%     legend('Green', 'Blue', 'Draw');

    hold on;
    plot(1:maxDraws, drawSequence(dri, :)*0.8 + 0.1, 'ro');

    axis([0 maxDraws + 1 0 1]);

end

figure;

R.sample = -0.005;

[r, Qsa] = backWardInduction(Ntrials, maxDraws, drawSequence, R);

for dri = 1 : 4
    subplot(2,2,dri);
    plot(1:maxDraws, squeeze(Qsa(dri, :, :)), 'LineWidth', 2);
%     legend('Green', 'Blue', 'Draw');

    choiceTrial = find(squeeze(Qsa(dri, :, 3)) - max(squeeze(Qsa(dri, :, 1:2))') < 0);
    text(choiceTrial(1), 0.1, 'X', 'horizontalalignment', 'center');


    hold on;
    plot(1:maxDraws, drawSequence(dri, :)*0.8 + 0.1, 'ro');

    axis([0 maxDraws + 1 0 1]);

end