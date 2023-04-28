% Function to generate a sequence of beads with the given ratio
function beads = generateBeadSequence(nBeads, blueRatio)

nBlue = round(nBeads * blueRatio);
nGreen = nBeads - nBlue;
beads = [ones(1,nBlue) ones(1,nGreen)*2];
beads = beads(randperm(nBeads));

end