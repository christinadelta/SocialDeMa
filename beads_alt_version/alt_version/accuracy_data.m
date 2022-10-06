function [all_accuracy, all_draws, all_points] = accuracy_data( session_data )

num_subs = 11;
if exist('block_types');
    num_conds = 4;
else
    num_conds = 1;
end;

iterator = 1;
for session=1:size( session_data(:,1) );
    
    clear block_type behavior blue_urn;
    load(session_data(session,:),'blue_urn','behavior','sequence');
    
    %2022
    %limit to blue urn or green urn trials
    blue_or_green = 0;  %one keeps blue urn and two keeps green urn - should match indices third dim of Q vector below (1 blue 2 green 3 sample) if meant to use to assess accuracy
    blue_urn_trials = find(blue_urn==blue_or_green);
    for i=1:numel(blue_urn_trials);
    
        temp_s{i} = sequence{blue_urn_trials(i)};
        temp_bu(i) = blue_or_green;
    
    end;
    
    allseq = sequence;
    blueurn = blue_urn;
    
    
    sequence = temp_s;
    blue_urn = temp_bu;
    
    num_trials = size(blue_urn,2);
    
    %     for trial = 1:num_trials;
    %for trial = 11:(num_trials-10);
    for trial = 1:num_trials;
        
        
        %data which is indexed by sequence/trial
        all_sessions_blue_urn( iterator ) = blue_urn( trial );
        try;    %in caes I didn't use different conditions
            all_sessions_block_types( iterator ) = block_type( trial );
        catch
            all_sessions_block_types( iterator ) = 1;   %just make them all 1
        end;
        all_sessions_behavior{ iterator } = behavior{trial};
        
        %need to find trial on which choice was made
        for option=1:size(behavior{trial},1);
            temp1 = behavior{trial};
            model_behavior{trial}(option,:) = temp1(option,:);  %need to kick out grey squares for modelling below
            if temp1(option,1) == 1 | temp1(option,2) == 1
                break
            end;
        end;
        all_sessions_draws( iterator ) = option - 1;
        
        %         all_sessions_draws( iterator ) = size(behavior{trial},1) - 1;
        iterator = iterator + 1;
    end;           %ends trials
    
end;  %ends sessions

for types = 1:num_conds;
    
    clear temp type_indices block_batches behavior_batches;
    type_indices = find(all_sessions_block_types == types);
    %temp=[];
    for example = 1:size( type_indices, 2);
        
        block_type_type(example) = all_sessions_block_types( type_indices(example) );
        blue_urn_type(example) = all_sessions_blue_urn( type_indices(example) );
        behavior_type{example,:} = all_sessions_behavior{ type_indices(example)}(:,1:3);    %all_sessions_behavior
        draws_type(example,:) = all_sessions_draws( type_indices(example) );
        %temp = [temp; behavior_type{example,:}];
    end;
    
    %initialize
    num_cell_trials = size( type_indices,2);
    correct = zeros( num_cell_trials, 1);
    
    for trial=1:num_cell_trials;
        
        if blue_urn_type(trial) == 1 && sum( behavior_type{trial}(:,1))>0;
            correct(trial,1) = 1;
        elseif blue_urn_type(trial) == 0 && sum( behavior_type{trial}(:,2))>0;
            correct(trial,1) = 1;
        end;
        
    end;    %loop thru trials this session
    
    all_accuracy(1) = mean(correct);
    all_draws(1) = mean( draws_type );
    all_points(1) = 40 + (sum(correct)*10) - sum(draws_type);
    
    %%%%%%%now evaluate ideal observer model matched to this subject
           for seq_num = 1: size(sequence, 2); %loop through all sequences
            if blue_urn(seq_num) == 0;
                seq_ones = find( sequence{seq_num} == 1);
                seq_twos = find( sequence{seq_num} == 2);
                sequence{seq_num}(seq_ones) = 2;
                sequence{seq_num}(seq_twos) = 1;
            end;    %ends test for green urn sequences
           end;    %ends loop through sequences
           
           %need sequence and behavior, but behavior needs to be in its old format
%            alpha = 1;  %softmax
%            Cw = 0;
%            q = 0.7;
%            Cs = 0;
%            [ll, pickTrial, dQvec, ddec, aQvec choice] = estimateLikelihoodf(alpha,Cw,q,Cs,cell2mat(sequence'),model_behavior,1);
%            
           R.correct = 10;
           R.error = 0;
           R.sample = 0;
           R.q = 0.7;
           
           seq_mat = cell2mat(sequence');
           seq_mat(find(seq_mat==2)) = 0;
           
           [r, Qsa] = backWardInduction(size(seq_mat,1), size(seq_mat,2), seq_mat, R);
           
           for dri=1:size(seq_mat,1);
               choiceTrial = find(squeeze(Qsa(dri, :, 3)) - max(squeeze(Qsa(dri, :, 1:2))') < 0);   %which options this seq have an urn > sample
               pickTrial(dri) = choiceTrial(1);    %which option was the first one where saample was inferior to urn (choice)?
               [ma ma_i] = max(squeeze(Qsa(dri,pickTrial(dri),:)));    %index of max value on choice option (which urn was chosen/best)?
               %2022
               if (ma_i == 1 & blue_or_green ==1) | (ma_i == 2 & blue_or_green == 0); 
                   choice(dri) = 1;
               else; choice(dri) = 0;
               end;
               
%                %2018
%                choice(dri) = ma_i;  %assign chosen urn
%                choice(find(choice==2)) = 0; %recode it so it can be summed
           end; %sequences loop
           
           all_accuracy(2) = mean(choice);
           all_draws(2) = mean(pickTrial);
           all_points(2) = 40 + (sum(choice)*10) - sum(pickTrial);
           
end; 

return