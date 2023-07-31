function modeloutput = fitAllModel_MC(R,thiscond_seqmat,cond_choices, urntype)

% run bayesbeads function 
% this function calls .m files located in model_fit folder 
% and stores model fit (and parameter recovery outputs

[minParams, ll, ~,~,~,~]    = bayesbeads_b(thiscond_seqmat, cond_choices, R);

modeloutput.fittedX         = minParams;
modeloutput.NLL             = ll;


end % end of function