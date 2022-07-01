function drawSequence = generateDrawSequences(q, maxDraws, Ntrials)

drawSequence = (rand(Ntrials, maxDraws) > q) + 1;

return