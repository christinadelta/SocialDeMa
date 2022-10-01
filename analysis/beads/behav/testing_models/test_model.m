% version 1 of model testing

% this script creates simulated data and runs ideal observer:
% 1. Nick's version
% 2. Bruno's version 

% both model appear to behave similarly when Cost_sample parameter is set to 0, Bruno's model draws more than
% Nick's model in the 0.6 probability condition when cost_sample is included!

clear all
clc

% define parameters 
R.correct = 10;
R.error = -10;
%R.Cw = -20;
R.sample =-0.25;
R.q = 0.6;
R.alpha = 1;

% generate data
for rs = 1 : 20
    if R.q == 0.8

        blue_b = ones(1,8);
        green_b = ones(1,2)*2;
    else
        blue_b = ones(1,6);
        green_b = ones(1,4)*2;
    end

    sq = cat(2,blue_b,green_b);

    tmp = randperm(length(sq));

    sequence(rs,:) = sq(tmp);
    
end

% create random 1s and 0s (for urn type)
temp_urns = cat(1, ones(10,1), zeros(10,1));

seq_mat = sequence;

% % if green urn switch indecies 
for s = 1:length(temp_urns)
    
    if temp_urns(s) == 0
        
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
    
    if (ma_i == 1 & temp_urns(dri) == 1) | (ma_i == 2 & temp_urns(dri) == 0)
        
        choice(dri) = 1;
    else
        choice(dri) = 0;
    end  
end %seque

all_acc = mean(choice);
all_picks = mean(pickTrial);
all_pts = (sum(choice==1)*R.correct) + (sum(choice == 0)*R.error) - (sum(pickTrial)*R.sample);
all_pts2 = (sum(choice==1)*R.correct) + (sum(choice == 0)*R.error) + (sum(pickTrial)*R.sample);

% now run nick's code (estimateLikelihood)
[ll, picktrl, dQvec, ddec, aQvec choices] = estimateLikelihoodf(seq_mat,R);

%choices(find(choices==2)) = 0;
for c = 1:Ntrials
    
    if (choices(c) == 1 & temp_urns(c) == 1) | (choices(c)== 2 & temp_urns(c) == 0)
        allchoices(c) = 1;
    else
        allchoices(c) = 0;
    end 
end

all_accuracy = mean(allchoices==1);
all_draws = mean(picktrl);
all_points = (sum(allchoices==1)*10) + (sum(allchoices==0)*-10) - (sum(picktrl)*-0.25);
all_points2 = (sum(allchoices==1)*10) + (sum(allchoices==0)*-10) + (sum(picktrl)*-0.25);

