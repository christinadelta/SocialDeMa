function p = PG(q, numDraws, numGreen)

p = 1/(1 + (q/(1-q))^(numDraws-2*numGreen));

end