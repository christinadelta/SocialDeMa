function [choiceStop, choiceCont, difVal, currentRnk] = computeSecretary(Generate_params, sampleSeries, prior, N, list, Cs,  minValue)

sdevs = 8;
dx = 2*sdevs*sqrt(prior.sig)/100;
x = ((prior.mu - sdevs*sqrt(prior.sig)) + dx : dx : ...
    (prior.mu + sdevs*sqrt(prior.sig)))';

Nchoices = length(list.vals);

Nconsider = length(sampleSeries);
if Nconsider > N
    Nconsider = N;
end

difVal      = zeros(1, Nconsider);
choiceCont  = zeros(1, Nconsider);
choiceStop  = zeros(1, Nconsider);
currentRnk  = zeros(1, Nconsider);





retrun