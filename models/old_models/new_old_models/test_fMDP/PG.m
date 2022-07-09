function p = PG(thisq, nd, ng)

p = 1/(1 + (thisq/(1-thisq))^(nd-2*ng));

return