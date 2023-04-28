function modeloutput = fitAllModel(R,thiscond_seqmat,cond_choices, urntype)

% run bayesbeads function 
% this function calls .m files located in model_fit folder 
% and stores model fit (and parameter recovery outputs

[minParams, ll, Qsad, cprob, model_samples,model_urnchoice]     = bayesbeads_b(thiscond_seqmat, cond_choices, R);

modeloutput.fittedX         = minParams;
modeloutput.NLL             = ll;
modeloutput.Qvals           = Qsad;
modeloutput.choiceProbs     = cprob;
modeloutput.modelSamples    = model_samples;

% compute average model samples
modeloutput.avSamples       = mean(model_samples);

% compute model performance 
for t = 1:size(model_urnchoice,1)
                
    if (model_urnchoice(t) == 1 & urntype(t) == 1) | (model_urnchoice(t) == 2 & urntype(t) == 0)
        model_choice(t) = 1;
    else
        model_choice(t) = 0;
    end % end of condition
end % end of trials loop

modeloutput.modelPerformance = mean(model_choice);

end % end of function