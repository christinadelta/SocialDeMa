% TESTING IDEAL OBSERVER MODELS WITH SIMULATED DATA

clear all
clc

% define parameters 
R.correct   = 10;
R.error     = -10;
%R.Cw       = -20;
R.sample    = -0.25;
R.alpha     = 1;
probs       = [0.8 0.6];

% generate data
for i = 1:2
    if i == 1
        for j = 1:20
            
            blue_b = ones(1,8);
            green_b = ones(1,2)*2;
            
            sq = cat(2,blue_b,green_b);
            tmp = randperm(length(sq));
            sequence_e(j,:) = sq(tmp);
 
        end % end of boots loop
        
    elseif i == 2
        
        for j = 1:20
            
            blue_b = ones(1,6);
            green_b = ones(1,4)*2;
            
            sq = cat(2,blue_b,green_b);
            tmp = randperm(length(sq));
            sequence_d(j,:) = sq(tmp);
        end % end of boots loop
    end % end of if statement
end % end of condition loop

% generate 20 random numbers from 1-10
drawse = [2 4 3 3 4 2 5 1 2 3 5 4 6 2 1 2 4 3 4 1];
drawsd = [4 3 5 5 7 4 5 3 5 5 6 4 4 5 7 5 4 4 3 4];

for cond = 1:2
    
    % create random 1s and 0s (for urn type)
    temp_urns = cat(1, ones(10,1), zeros(10,1));
    rand_urns = randperm(length(temp_urns));
    urns(:,cond) = temp_urns(rand_urns);
    R.q = probs(cond);
    
    % which sequences to use? which condition are we in?
    if cond == 1
        seq_mat = sequence_e;
        condraws = drawse;
    else
        seq_mat = sequence_d;
        condraws = drawsd;
    end

    % create choice vectors for each condition/trial
    for i = 1:length(urns)
        thisdraw = condraws(1,i);
        
        for tr = 1:thisdraw
            if tr ~= thisdraw % if this is not the last draw
                t(tr,1:2) = 0;
                t(tr,3) = 1;
            elseif urns(i) == 1
                t(tr,1) = 1
                t(tr,2:3) = 0;
            elseif urns(i) == 0
                t(tr,1) = 0
                t(tr,2) = 1
                t(tr,3) = 0
            end
        end
        
        choiceVec{cond,i} = t;
        clear t
    end % end of choice vec 
    
    % % if green urn switch indecies 
    for s = 1:length(urns)
        
        if urns(s) == 0
            
            seq_ones = find( seq_mat(s,:) == 1);
            seq_twos = find( seq_mat(s,:) == 2);
            seq_mat(s,seq_ones) = 2;
            seq_mat(s,seq_twos) = 1;
            
        end
    end
    
    % this is for the model
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
        
        if (ma_i == 1 & urns(dri) == 1) | (ma_i == 2 & urns(dri) == 0)
            
            choice(dri) = 1;
        else
            choice(dri) = 0;
        end
    end %sequ
    
    all_acc(1,cond) = mean(choice)
    all_picks(1,cond) = mean(pickTrial)
    all_pts(1,cond) = (sum(choice==1)*R.correct) + (sum(choice == 0)*R.error) + (sum(pickTrial)*R.sample)
    
    % now run nick's code (estimateLikelihood)
    [ll, picktrl, dQvec, ddec, aQvec choices] = estimateLikelihoodf(seq_mat,R);
    
    %choices(find(choices==2)) = 0;
    for c = 1:Ntrials
        
        if (choices(c) == 1 & urns(c) == 1) | (choices(c)== 2 & urns(c) == 0)
            allchoices(c) = 1;
        else
            allchoices(c) = 0;
        end
    end
    
    all_accuracy(1,cond) = mean(allchoices==1)
    all_draws(1,cond) = mean(picktrl)
    all_points(1,cond) = (sum(allchoices==1)*10) + (sum(allchoices==0)*-10) + (sum(picktrl)*-0.25)

 
end % emd of condition loop

