function p = PG_io(q, numDraws, numGreen)

% "compute conditional probabilities of the majority urn" (green? blue?)
% This function runs at each draw and needs: number of draws so far, number of green beads (or
% blue beads) and the probability of drawing a majority colour bead (0.8,
% 0.6) 

% p(g|nd,ng):

p = 1/(1 + (q/(1-q))^(numDraws-2*numGreen));

end