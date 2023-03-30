function p = PG(q, nd, ng)

p = 1/(1 + (q/(1-q))^(nd-2*ng));

return