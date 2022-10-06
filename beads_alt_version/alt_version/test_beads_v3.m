% playing around with beads stuff


condition_dirs = {'/Users/christinadelta/Desktop/ideal_observer/alt_version/combined/domain01','/Users/christinadelta/Desktop/ideal_observer/alt_version/combined/domain02'};
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