% This script runs ideal observer and model fitting

% adds the ideal observer and model fitting dirs
currentpath     = pwd;
biobserverpath  = fullfile(pwd, 'brunos_io');
iobserverpath   = fullfile(pwd, 'ideal_observer');

bmodelfitpath   = fullfile(pwd, 'brunos_modelfit');
modelfitpath    = fullfile(pwd, 'model_fitting');



clear all 
clc
% load workspace 
load('workspace3.mat')

%% run ideal observer

% first add modelpath to the path
addpath(genpath(biobserverpath)); 
addpath(genpath(iobserverpath)); 

% define parameters of the ideal observer
R.alpha             = 1;            % softmax stochasticity parameter (for fitting to human behaviour) - this is not needed here
R.error             = -10;          % cost for being wrong
%     R.diff              = -20;          % The difference between the rewards for being correct (in this case no reward 10) and the cost of being wrong (-10).
R.correct           = 10;           % reward for being correct
R.q                 = [0.8 0.6];    % proportion of the majority value in sequence (60:40 split in this case)
R.sample            = -0.25;            % the cost to sample

for sub = 1:nsubs

    for cond = 1:conditions
        
        
        thiscond_data    = cond_data{1,sub}{1,cond};
        thiscond_seq     = allsub_sequences{1,sub}{1,cond};
        
        % extract sequences from cell and store in matrix (26x10)
        for i = 1:size(thiscond_seq,2)
            thiscond_seqmat(i,:) = thiscond_seq{1,i};
        
        end
        
        tmp_seqmat = thiscond_seqmat;
        
        % what is the probability of this cond? 
        if cond == 1
            thisq = R.q(1);
        else 
            thisq = R.q(2);
        end
        
        R.thisq = thisq;
        
        %%%% First run Bruno's io 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % what is the urntype?
        urntype = thiscond_data(:,4);
        
        for u = 1:length(urntype)
            if urntype(u) == 0 % if green urn switch index coding
                seq_ones = find(thiscond_seqmat(u,:) == 1);
                seq_twos = find(thiscond_seqmat(u,:) == 2);
                thiscond_seqmat(u,seq_ones) = 2;
                thiscond_seqmat(u,seq_twos) = 1;
            end 
        end
        
        % recode 2s to 0s for backward induction 
        thiscond_seqmat(find(thiscond_seqmat==2))=0;
        
        % run backward induction (bruno's code)
        [r, Qsat] = backWardInduction(thiscond_seqmat, R);
        
        % store ideal observer output
        bio_output(cond).r       = r;
        bio_output(cond).Qsat    = Qsat;
        
        % loop over condition trials to compute choices, picktrials and acc
        for i = 1: totaltrials/2
            
            choiceTrial             = find(squeeze(Qsat(i,:,3)) - max(squeeze(Qsat(i,:,1:2))') < 0); % which options this trial vec have an urn > sample
            pickTrial(i)            = choiceTrial(1); % pick the first of the choices
            [ma ma_i]               = max(squeeze(Qsat(i,pickTrial(i),:))); % which of the two urn was chosen? [based on AQ values]
            
            if (ma_i == 1 & urntype(i) == 1) | (ma_i == 2 & urntype(i) == 0)
                choice(i) = 1;
            else
                choice(i) = 0;
            end
        end  
        
        % for each subject model instance and condition, calculate acc and
        % draws
        allsub_bioacc(sub,cond)     = mean(choice);
        allsubs_biodraws(sub,cond)  = mean(pickTrial); 
        allsubs_biopoints(sub,cond) = (sum(choice==1)*R.correct) + (sum(choice==0)*R.error) + (sum(pickTrial)*R.sample);
        
        
        %%%%%%%%%%%%%  run Nick's io version
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % run estimateLikelihoodf (Nick's code)
        [ll, picktrl, dQvec, ddec, aQvec choices] = estimateLikelihoodf_io(tmp_seqmat,R);
        
        io_output(cond).aQvec       = aQvec;
        io_output(cond).picktrl     = picktrl;
        
        choices(find(choices==2))   = 0; % re-code incorrect responses
        
        % store all io outpus
        allsub_ioacc(sub,cond)     = mean(choices==1);
        allsubs_iodraws(sub,cond)  = mean(picktrl); 
        allsubs_iopoints(sub,cond) = (sum(choices==1)*R.correct) + (sum(choices==0)*R.error) + (sum(picktrl)*R.sample);
 
  
    end % end of condition loop
    
    clear r Qsat thiscond_data thiscond_seq thiscond_seqmat choice pickTrial urntype
 
end % end of subs loop

clear R 

%% Run model fitting 

% % IN CASE YOU WANT TO CHECK MODEL FITTING TOO
% 
% % first add model fitting to the path
% addpath(genpath(bmodelfitpath));
% addpath(genpath(modelfitpath));
%     
% 
% for sub = 1:nsubs
%     
%     % define parameters of the ideal observer
%     R.initbeta          = 1;            % softmax stochasticity parameter (for fitting to human behaviour) - this is not needed here
%     R.error             = -10;          % cost for being wrong
%     R.diff              = -20;          % The difference between the rewards for being correct (in this case no reward 10) and the cost of being wrong (-10).
%     R.correct           = 10;           % reward for being correct
%     R.q                 = [0.8 0.6];    % proportion of the majority value in sequence (60:40 split in this case)
%     R.initsample        = -0.25;        % the cost to sample
%     
%     for cond = 1:conditions
%         
%         thiscond_data       = cond_data{1,subI}{1,cond};
%         thiscond_seq        = allsub_sequences{1,sub}{1,cond};
%         thiscond_choiceVec  = allsub_choiceVec{1,sub}{1,cond};
%         
%         % what is the probability of this cond? 
%         if cond == 1
%             thisq = R.q(1);
%         else 
%             thisq = R.q(2);
%         end
%         
%         R.thisq = thisq;
%         
%         % extract sequences from cell and store in matrix (26x10)
%         for i = 1:size(thiscond_seq,2)
%             thiscond_seqmat(i,:) = thiscond_seq{1,i};
%         end
%         
%         % what is the urntype?
%         urntype = thiscond_data(:,4);
% 
%         for u = 1:length(urntype)
%             if urntype(u) == 0 % if green urn switch index coding
%                 seq_ones = find(thiscond_seqmat(u,:) == 1);
%                 seq_twos = find(thiscond_seqmat(u,:) == 2);
%                 thiscond_seqmat(u,seq_ones) = 2;
%                 thiscond_seqmat(u,seq_twos) = 1;
%             end 
%         end
% 
%         % recode 2s to 0s for backward induction 
%         thiscond_seqmat(find(thiscond_seqmat==2))=0;
%         
%         %%%%%%%%%%%%%% fit Nick's version of the model
%         [mparams, lla, aQvec]           = bayesbeads(thiscond_seqmat, thiscond_choiceVec, R);
%         
%         % store model-fitting output
%         allsubs_cs_mn(subI,cond)        = mparams(1);
%         allsubs_beta_mn(subI,cond)      = mparams(2);
%         allsubs_lla_mn(subI,cond)       = lla;
%         allsubs_AQs_mn{1,subI}          = aQvec;
%         
%         
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         
%         %%%%%%% fit Bruno's version of the model
%         [minParams, ll, Qsad, cprob]    = bayesbeads_b(thiscond_seqmat, thiscond_choiceVec, R);
%         
%         allsubs_cs_mb(subI,cond)        = minParams(1);
%         allsubs_beta_mb(subI,cond)      = minParams(2);
%         allsubs_ll_mb(subI,cond)        = ll;
%         allsubs_AQs_mb{1,subI}          = Qsad;
%         
%   
%     end % end of conditions loop
%     
% end % end of subs loop

