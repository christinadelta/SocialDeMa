% TESTING IDEAL OBSERVER MODELS AND MODEL-FITTING WITH SIMULATED DATA

% this script simulates data for N participants and runs:
% 1) Ideal observer (Nick's and Bruno's version)
% 2) Model fitting (Nick's and Bruno's version)

% THE MAIN PURPOSE FOR THIS IS:
% 1. TO GO THROUGH THE MODEL FITTING PROCEDURE IN DETAIL
% 2. TO GO THROUGH NICK'S "ESTIMATELIKELIHOOD" APPROACH AND LOOK WHY THERE
% ARE DIFFERENCES WITH BRUNO'S "BACKWARDSINDUCTION" APPROACH 

% Data to be simulated:
% 1) sequences (in 2 conditions)
% 2) subject draws 
% 3) choice vectors 

% Created: 19/10/2022

clear all
clc

startpath   = pwd;
modelfitb   = fullfile(startpath, 'model_fitting_bio');
addpath(genpath(modelfitb));

%% define initial vars 

% define parameters 
R.correct   = 10;
R.error     = -10;
%R.Cw       = -20;
R.sample    = -0.25;
R.alpha     = 1;
probs       = [0.8 0.6];
nsubs       = 3;

%% loop over subjects and start simulatung data

for sub = 1:nsubs 
    
    % 1. generate sequence data
    for i = 1:2 % conditions
        if i == 1
            for j = 1:20 % number of sequences per condition

                blue_b = ones(1,8);
                green_b = ones(1,2)*2;

                sq = cat(2,blue_b,green_b);
                tmp = randperm(length(sq));
                sequence_e(j,:) = sq(tmp);
                
                % add all sequences in 1 cell
                subseq{1, i} = sequence_e;

            end % end of j loop

        elseif i == 2

            for j = 1:20

                blue_b = ones(1,6);
                green_b = ones(1,4)*2;

                sq = cat(2,blue_b,green_b);
                tmp = randperm(length(sq));
                sequence_d(j,:) = sq(tmp);
                subseq{1, i} = sequence_d;
                
            end % end of boots loop
        end % end of if statement 
    end % end of condition loop
    
    clear i j sequence_d sequence_e tmp sq blue_b green_b
    
    % 2. for each condition (and sequenc) generate 20 random numbers from 1-10
    for i = 1:2
        
        xmin            = 1;
        xmax            = 10;
        n               = 20;
        x               = round(xmin + rand(1,n) * (xmax - xmin))';
        
        alldraws(:,i)   = x; % the draws include urn choice response (this will be subtracted later)
        
        clear xmin xmax x 
    end % end of cond loop
    
    % 3. create urntypes and vectors with choices 
    for i = 1:2
        
        % create random 1s and 0s (for urn type)
        temp_urns = cat(1, ones(10,1), zeros(10,1));
        rand_urns = randperm(length(temp_urns));
        urns(:,i) = temp_urns(rand_urns);
        
        clear temp_urns rand_urns
    
        % create choice vectors for each condition/trial
        for j = 1:length(urns)
            thisdraw = alldraws(j,i);

            for tr = 1:thisdraw
                if tr ~= thisdraw % if this is not the last draw
                    t(tr,1:2) = 0;
                    t(tr,3) = 1;
                elseif urns(j,i) == 1
                    t(tr,1) = 1
                    t(tr,2:3) = 0;
                elseif urns(j,i) == 0
                    t(tr,1) = 0
                    t(tr,2) = 1
                    t(tr,3) = 0
                end
            end

            choiceVec{i,j} = t;
            clear t
        end % end of choice vec 
    end % end of cond loop
    
    clear i j n thisdraw tr
    
    %% run ideal observer models 
    
    for i = 1:2
        
        % which sequences to use? which condition are we in?
        seq_mat     = subseq{1,i};
        condraws    = alldraws(:,i);
        condurns    = urns(:,i);
        
        cond_seq    = seq_mat; % this will be used in Nick's ideal observer model
        
        R.q     = probs(i);
        
        % % if green urn switch indecies 
        for s = 1:length(condurns)

            if condurns(s) == 0

                seq_ones = find( seq_mat(s,:) == 1);
                seq_twos = find( seq_mat(s,:) == 2);
                seq_mat(s,seq_ones) = 2;
                seq_mat(s,seq_twos) = 1;

            end
        end

        % this is for the model (the way that the majority colour beads are
        % calculated. If left as 1s and 2s then, then the number of
        % majority beads is misscalculated
        seq_mat(find(seq_mat==2)) = 0;
        
        Ntrials = size(seq_mat,1);
        maxDraws = size(seq_mat,2);

        % run backward induction (Bruno's code)
        [r, Qsat] = backWardInduction(Ntrials, maxDraws, seq_mat, R);

        for dri = 1:Ntrials

            choiceTrial = find(squeeze(Qsat(dri, :, 3)) - max(squeeze(Qsat(dri, :, 1:2))') < 0);   %which options this seq have an urn > sample
            pickTrial(dri) = choiceTrial(1);    %which option was the first one where saample was inferior to urn (choice)?
            [ma ma_i] = max(squeeze(Qsat(dri,pickTrial(dri),:)));    %index of max value on choice option (which urn was chosen/best)?
            %     choice(dri) = ma_i;  % assign chosen urn
            % choice(find(choice==2)) = 0; %recode it so it can be summed

            if (ma_i == 1 & condurns(dri) == 1) | (ma_i == 2 & condurns(dri) == 0)

                choice(dri) = 1;
            else
                choice(dri) = 0;
            end
        end %sequ
        
        % get averages and model points 
        all_acc_bio(sub,i)      = mean(choice);
        all_picks_bio(sub,i)    = mean(pickTrial);
        all_pts_bio(sub,i)      = (sum(choice==1)*R.correct) + (sum(choice == 0)*R.error) + (sum(pickTrial)*R.sample);
        
        % now run nick's code (estimateLikelihood)
        [ll, picktrl, dQvec, ddec, aQvec choices] = estimateLikelihoodf(cond_seq,R);
        
        choice(find(choice==2)) = 0; % re-code incorrect responses
        
        % get averages
        all_picks_io(sub,i)     = mean(picktrl);
        all_acc_io(sub,i)       = mean(choice==1);
        all_points_io(sub,i)    = (sum(choice==1)*10) + (sum(choice==0)*-10) + (sum(picktrl)*-0.25);
        
        clear ll picktrl dQvec ddec aQvec choices ma ma_i r Qsat dri correct choice k condraws condurns cond_seq seq_mat seq_ones 
        clear seq_twos Ntrials maxDraws pickTrial choiceTrial s 

    end % end of conditions loop
    
    clear i 
    
    %% run model fitting %% 
    
    % for now run only on Bruno's model version as it seems to work better
    % than Nick's 
    for i = 1:2
        
        % update R strcut
        R.q = probs(i);
        R.params = [-0.25; 1]; % [cost-to-sample and beta value]
        
        % which sequences to use? which condition are we in?
        seq_mat     = subseq{1,i};
        condraws    = alldraws(:,i);
        condurns    = urns(:,i);

        cond_seq    = seq_mat; % this will be used in Nick's ideal observer model
        
        R.q     = probs(i);
        
        % % if green urn switch indecies 
        for s = 1:length(condurns)

            if condurns(s) == 0

                seq_ones = find( seq_mat(s,:) == 1);
                seq_twos = find( seq_mat(s,:) == 2);
                seq_mat(s,seq_ones) = 2;
                seq_mat(s,seq_twos) = 1;

            end
        end

        % this is for the model (the way that the majority colour beads are
        % calculated. If left as 1s and 2s then, then the number of
        % majority beads is misscalculated
        seq_mat(find(seq_mat==2)) = 0;
        
        R.cond = i;
        
        % fit model 
        [minParams,ll,Qsad,cprob] = fmdpBeadsTest(seq_mat,choiceVec,condraws,R);
        
        % store results
        all_betas(sub,i) = minParams(2);
        all_cs(sub,i) = minParams(1);
        all_ll(sub,i) = ll;

    end
 
end % end of subjects loop












