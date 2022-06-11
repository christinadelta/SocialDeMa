%%% MODEL FITTING BEADS TASK

% sort behaviour

% run ideal observer 
[ll, pickTrial, dQvec, ddec, aQvec choice] = estimateLikelihoodf(alpha,Cw,thisq,Cs,this_sequence,1);

pick_trials{sub, block, trl}    = pickTrial;
blockdQvec{sub, block, trl}     = dQvec;
bloxkaQvec{sub, block, trl}     = aQvec;
blockchoices(sub, block, trl)   = choice;
blockddec{sub, block, trl}      = ddec;



