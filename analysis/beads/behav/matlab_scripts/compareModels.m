function [group_best best] = compareModels(R, cond_data,nsubs,allsub_choiceVec,allsub_sequences)

% Fits the models, computes BIC and select model based on BIC value
% To Compute the BIC value: Use the formula for BIC to calculate the BIC value for each model:


%            BIC = -2 * total NLL + k * log(n)   ---- at the group level
%            BIC = -2 * log-likelihood + k * log(n) --- at the subject level

% Where:
% 1. log-likelihood: The logarithm of the likelihood for the model given the data.
% 2. k: The number of parameters in the model.
% 3. n: The number of data points in the dataset (sequences).
% 4. total NLL: the su of NLL across participants 

% NOTE: 
% at the group and subjects levels n is not the same:
% n = number of subjects (at group level)
% n = number of sequences per subject (at the subject level)

% Now by datapoints... should it be number of draws per sequence or number
% of sequences? -- Number of sequences (asked chatGPT)

% MODEL COMPARISON WILL BE COMPLETED AT:
% 1. THE GROUP LEVEL 
% 2. SUBJECT LEVEL 

% -------------------------
%% fit model 1 

model           = 1; % which model is it?
R.freeparams    = 1; % how many free parameters?
totalNLL_easy   = 0;
totalNLL_diff   = 0;


% define free parameters
R.initsample    = R.Cs; % this won't be used in this model 
R.initbeta      = R.beta; 
R.model         = model;

% loop over subjects 
for sub = 1:nsubs

    subdata             = cond_data{1,sub};         % all data matrix
    sub_choicesvec      = allsub_choiceVec{1,sub};  % subject choices (two 26 by 3 matricies)
    sub_sequence        = allsub_sequences{1,sub};  % sequence that this-subject was presented with

    % loop over conditions 
    for cond = 1:2

        % what is the probability of this condition?
        if cond == 1

            thisq       = R.q(1);
        else 
            thisq       = R.q(2);
        end
    
        R.thisq         = thisq;

        % extract condition data
        sub_cond        = subdata{1,cond};
        cond_choices    = sub_choicesvec{1,cond};
        cond_sequence   = sub_sequence{1,cond};

        % extract sequences from cell and store in matrix (26x10)
        for s = 1:size(cond_sequence,2)
            thiscond_seqmat(s,:)                        = cond_sequence{1,s};
        end
        
        % what is the urntype?
        urntype                                         = sub_cond(:,4);

        for u = 1:length(urntype)
            if urntype(u) == 0 % if green urn switch index coding
                seq_ones                                = find(thiscond_seqmat(u,:) == 1);
                seq_twos                                = find(thiscond_seqmat(u,:) == 2);
                thiscond_seqmat(u,seq_ones)             = 2;
                thiscond_seqmat(u,seq_twos)             = 1;
            end 
        end

        % recode 2s to 0s for backward induction 
        thiscond_seqmat(find(thiscond_seqmat==2))       = 0;

        modeloutput                         = fitAllModel_MC(R,thiscond_seqmat,cond_choices,urntype);

        betaModelNLL(sub,cond)              = modeloutput.NLL;
        betamodelXval(sub,cond)             = modeloutput.fittedX;
        
        % sum NLL for all participants 
        if cond == 1
            totalNLL_easy                   = totalNLL_easy + modeloutput.NLL; % sum NLL for all participants
        elseif cond == 2
            totalNLL_diff                   = totalNLL_diff + modeloutput.NLL; % sum NLL for all participants
        end


    end % end of conditions loop
end % end of subjects loop

%% fit model 2

model           = 2; % which model is it?
R.freeparams    = 2; % how many free parameters?
totalNLL2_easy   = 0;
totalNLL2_diff   = 0;


% define free parameters
R.initsample    = R.Cs; % this won't be used in this model 
R.initbeta      = R.beta; 
R.model         = model;

% loop over subjects 
for sub = 1:nsubs

    subdata             = cond_data{1,sub};         % all data matrix
    sub_choicesvec      = allsub_choiceVec{1,sub};  % subject choices (two 26 by 3 matricies)
    sub_sequence        = allsub_sequences{1,sub};  % sequence that this-subject was presented with

    % loop over conditions 
    for cond = 1:2

        % what is the probability of this condition?
        if cond == 1

            thisq       = R.q(1);
        else 
            thisq       = R.q(2);
        end
    
        R.thisq         = thisq;

        % extract condition data
        sub_cond        = subdata{1,cond};
        cond_choices    = sub_choicesvec{1,cond};
        cond_sequence   = sub_sequence{1,cond};

        % extract sequences from cell and store in matrix (26x10)
        for s = 1:size(cond_sequence,2)
            thiscond_seqmat(s,:)                        = cond_sequence{1,s};
        end
        
        % what is the urntype?
        urntype                                         = sub_cond(:,4);

        for u = 1:length(urntype)
            if urntype(u) == 0 % if green urn switch index coding
                seq_ones                                = find(thiscond_seqmat(u,:) == 1);
                seq_twos                                = find(thiscond_seqmat(u,:) == 2);
                thiscond_seqmat(u,seq_ones)             = 2;
                thiscond_seqmat(u,seq_twos)             = 1;
            end 
        end

        % recode 2s to 0s for backward induction 
        thiscond_seqmat(find(thiscond_seqmat==2))       = 0;

        modeloutput                         = fitAllModel_MC(R,thiscond_seqmat,cond_choices,urntype);

        betaCsModelNLL(sub,cond)            = modeloutput.NLL;
        betaCsmodelX1val(sub,cond)          = modeloutput.fittedX(1);
        betaCsmodelX2val(sub,cond)          = modeloutput.fittedX(2);

        % sum NLL for all participants 
        if cond == 1
            totalNLL2_easy                   = totalNLL2_easy + modeloutput.NLL; % sum NLL for all participants
        elseif cond == 2
            totalNLL2_diff                   = totalNLL2_diff + modeloutput.NLL; % sum NLL for all participants
        end


    end % end of conditions loop
end % end of subjects loop

%% compute BIC values at group level 

% for each model, compute BIC for each condition 

k               = 1; % number of free params
k2              = 2;
n               = nsubs; 

groupBIC(1)     = -2 * totalNLL_easy + k * log(n);      % beta model easy 
groupBIC(2)     = -2 * totalNLL_diff + k * log(n);      % beta model difficult
groupBIC(3)     = -2 * totalNLL2_easy + k2 * log(n);    % beta and Cs model easy 
groupBIC(4)     = -2 * totalNLL2_diff + k2 * log(n);    % beta and Cs model difficult


%% Compute BIC values at subject level 

% now lets compute BIC values for each model and condition  for each
% subject seperately 
k               = 1; % number of free params
k2              = 2;
n               = 26; % number of sequences 

% loop over subjects 
for sub = 1:nsubs

    % compute BIC vals for each model and condition
    BIC(sub,1) = -2 * betaModelNLL(sub,1) + k * log(n);     % beta model easy 
    BIC(sub,2) = -2 * betaModelNLL(sub,2) + k * log(n);     % beta model difficult
    BIC(sub,3) = -2 * betaCsModelNLL(sub,1) + k2 * log(n);  % beta and Cs model easy 
    BIC(sub,4) = -2 * betaCsModelNLL(sub,2) + k2 * log(n);  % beta and Cs model difficult

end % end of subjects loop

%% Find best model at group level

[m, ibest]          = min(groupBIC);
group_best          = ibest;

% best = groupBIC == m;
% best        = best / sum(best);

%% Find best model at subject level

% loop over subjects 
for sub = 1:nsubs

    tempBIC             = BIC(sub,:);
    [m, ibest]          = min(tempBIC);
    best(sub,1)         = ibest;

end % end of subjects loop

end % end of function 