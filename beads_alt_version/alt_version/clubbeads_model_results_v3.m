function [] = analyze_results;

%v3: playing around with this code in 2022. v2 was for club beads in 2018

disp('analyzing behavioral data');
%
% cd E:\beads\old_cluster_stuff\fmri;
% cd C:\matlab_files\fiance;

condition_dirs = {'C:\matlab_files\fiance\club_beads\data_171218\combined\domain01','C:\matlab_files\fiance\club_beads\data_171218\combined\domain02'};
num_conds = size(condition_dirs,2);

% temp_data = spm_select(Inf, 'mat', 'select subjects');
% num_subs = size(temp_data,1);

graph_font = 14;
% num_subs = 17;
% num_conds = 4;
% all_draws = cell(num_subs,1);
% all_accuracy = cell(num_subs,1);
all_draws = [];
all_accuracy = [];
% all_session_draws = cell(num_subs,1);
% ave_draws = zeros(num_subs, num_conds);
% ave_accuracy = zeros(num_subs, num_conds);
% fig_name_draws = 'draws_graph';
% fig_name_accuracy = 'accuracy_graph';
% fig_name_drawmatrix = 'draws_matrix';
% fig_name_accuracymatrix = 'accuracy_matrix';


iterator = 1;
for this_condition = 1:num_conds;
    
    subs = dir([condition_dirs{this_condition} filesep 'phase2_sequence_data_clubbeads_sub*_domain*.mat']);
    num_subs = size(subs,1);
    
    for subject = 1:num_subs;
        
        session_data = [condition_dirs{this_condition} filesep subs(subject).name];

        %v2: now extracts accuracy and draws from raw data using this funct
        [all_accuracy( this_condition,subject,:), all_draws( this_condition,subject,: ) all_points( this_condition,subject,: )] ...
            = accuracy_data( session_data );

    end;    %subjects
end;    %conditions

draws = squeeze(nanmean(all_draws,2));
correct = squeeze(nanmean(all_accuracy,2));
points = squeeze(nanmean(all_points,2));

draws_ci = (squeeze(nanstd(all_draws,1,2))./sqrt(size(all_draws,2))).*1.96;
correct_ci = (squeeze(nanstd(all_accuracy,1,2))./sqrt(size(all_accuracy,2))).*1.96;
points_ci = (squeeze(nanstd(all_points,1,2))./sqrt(size(all_points,2))).*1.96;

figure; 
subplot(1,3,1);
% barweb( draws', draws_ci');
bar( draws');
axis square; box off;
set(gca, 'xticklabel', {'Participants', 'Model'}, 'FontName', 'Ariel', 'FontSize', graph_font, 'FontWeight', 'bold') ;
xlabel('Choice domain'); ylabel('number draw choices');
set(gcf, 'Color', [1 1 1]); ylim([0 10]); legend('Club','Classic'); legend boxoff;

subplot(1,3,2);
% barweb( correct', correct_ci');
bar(correct');
axis square; box off;
set(gca,'ytick',[0:.10:1.00], 'xticklabel', {'Participants', 'Model'}, 'FontName', 'Ariel', 'FontSize', graph_font, 'FontWeight', 'bold') ;
xlabel('Choice domain'); ylabel('proportion correct choices');
set(gcf, 'Color', [1 1 1]); ylim([0 1.1]); 

subplot(1,3,3);
% barweb( points', points_ci');
bar( points');

axis square; box off;
set(gca, 'xticklabel', {'Participants', 'Model'}, 'FontName', 'Ariel', 'FontSize', graph_font, 'FontWeight', 'bold') ;
xlabel('Choice domain'); ylabel('points');
set(gcf, 'Color', [1 1 1]); 
%ylim([0 1.1]); 
% saveas( gcf, fig_name_draws, 'fig' );

fprintf(' ');

%only works if equal Ns!!! Careful it will create fae subs full of zeros if
%one cond has < N

%create an output file that can be imported into SPSS - Also the graphs
%above may be incorrect too!
club_human_draws = squeeze(all_draws(1,:,1))';
classic_human_draws = squeeze(all_draws(2,:,1))';
club_model_draws = squeeze(all_draws(1,:,2)');
classic_model_draws = squeeze(all_draws(2,:,2))';
human_draws = [[ones(size(club_human_draws,1),1); zeros(size(classic_human_draws,1),1)] [club_human_draws; classic_human_draws]];
model_draws = [club_model_draws; classic_model_draws];

club_human_accuracy = squeeze(all_accuracy(1,:,1))';
classic_human_accuracy = squeeze(all_accuracy(2,:,1))';
club_model_accuracy = squeeze(all_accuracy(1,:,2)');
classic_model_accuracy = squeeze(all_accuracy(2,:,2))';
human_accuracy = [club_human_accuracy; classic_human_accuracy];
model_accuracy = [club_model_accuracy; classic_model_accuracy];

club_human_points = squeeze(all_points(1,:,1))';
classic_human_points = squeeze(all_points(2,:,1))';
club_model_points = squeeze(all_points(1,:,2)');
classic_model_points = squeeze(all_points(2,:,2))';
human_points = [club_human_points; classic_human_points];
model_points = [club_model_points; classic_model_points];


out_matrix = [human_draws model_draws human_accuracy model_accuracy human_points model_points];
dlmwrite('project2019_club_beads.txt',out_matrix,' ');

disp('audi5000')


%%
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
        
        if blue_urn_type(trial) == 1 & sum( behavior_type{trial}(:,1))>0;
            correct(trial,1) = 1;
        elseif blue_urn_type(trial) == 0 & sum( behavior_type{trial}(:,2))>0;
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
           R.sample = -1;
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
           
end;    %ends types

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% for running synthetic daata
function [reward, Qsat] = backWardInduction(Ntrials, maxDraws, drawSequence, R)

K = 3;

reward = zeros(Ntrials, 1);

Qsat = zeros(Ntrials, maxDraws, K);

parfor trial = 1 : Ntrials
    
    Qsad = zeros(maxDraws, 3);
    
    for draw = 1 : maxDraws
                   
        Qsad(draw, :) = backWardUtility(drawSequence(trial, :), draw, maxDraws, R);
        
    end
    
    Qsat(trial, :, :) = Qsad;
    
    %%% randomize choice for symmetric values
    Qsac = Qsad + 0.000001*randn(maxDraws, K);
    
    Qsa1 = Qsac(:, 1) - Qsac(:, 3);
    Qsa2 = Qsac(:, 2) - Qsac(:, 3);
    
    choice1 = find(Qsa1 > 0);
    choice2 = find(Qsa2 > 0);
    
    if isempty(choice1)
        choice1(maxDraws+1) = 1;
    end
    
    if isempty(choice2)
        choice2(maxDraws+1) = 1;
    end
    
    if choice1(1) < choice2(1)
        reward(trial) = 1;
    else
        reward(trial) = 0;
    end
    
            
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Qsa = backWardUtility(drawSequence, draw, maxDraws, R)

utility = zeros(maxDraws, maxDraws+1);

ng = sum(drawSequence(1:draw));

for drawi = maxDraws : -1 : (draw + 1)
        
    [utility] = stateUtilityBeads(utility, drawi, draw, maxDraws, ng, R);
    
end
    
Qsa = actionValueBeads(utility, R, draw, ng, draw, maxDraws);
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function utility_t = stateUtilityBeads(utility, drawi, draw, maxDraws, ng, R)

utility_t = zeros(maxDraws, maxDraws+1);

futureDraws = drawi - draw;

ndf = drawi;

for greenDraws = 0 : futureDraws
    
    ngf = ng + greenDraws;

    Qsa = actionValueBeads(utility, R, ndf, ngf, drawi, maxDraws);

    utility_t(ndf, ngf+1) = max(Qsa);        
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Qsa = actionValueBeads(utility, R, nd, ng, drawi, maxDraws)

pg = PG(R.q, nd, ng);

pb = 1 - pg;

QG = R.correct*pg + R.error*pb;
QB = R.correct*pb + R.error*pg;

if drawi < maxDraws

    QD = R.sample + pb*((1-R.q)*utility(nd+1, ng+1+1) +   (R.q)*(utility(nd+1, ng+1))) + ...
                    pg*(  (R.q)*utility(nd+1, ng+1+1) + (1-R.q)*(utility(nd+1, ng+1)));
                
else
    
    QD = 0;
    
end

Qsa = [QG; QB; QD];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p = PG(q, nd, ng)

p = 1/(1 + (q/(1-q))^(nd-2*ng));
