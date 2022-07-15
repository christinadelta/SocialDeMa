function [num_samples ranks choiceStop_all choiceCont_all] = generate_a_models_data(Generate_params)

% returns subject*sequences matrices of numbers of draws and ranks

% note: What is called which_model here is param_to_fit(model) outside this programme,
% at the base level

this_sub = 1; % Need to assign each sub to output array by how many have been run, rathere than by sub num
for num_subs_found = Generate_params.num_subs_to_run
    
    if numel(Generate_params.num_subs_to_run) > 1 %i.e., if model fitting to a single subject is not going on here
        disp(...
            sprintf('generating performance for preconfigured modeli %d name %s subject %d' ...
            , Generate_params.current_model ...
            ,Generate_params.model( Generate_params.current_model ).name ...
            , num_subs_found ...
            ) );
    end
    
    for sequence = 1:Generate_params.num_seqs
        
        %         if Generate_params.model(Generate_params.current_model).log_or_not == 1;
        %             Generate_params.binEdges_reward = ...
        %                 linspace(log(Generate_params.rating_bounds(1)),log(Generate_params.rating_bounds(2)),Generate_params.nbins_reward+1);   %organise bins by min and max
        %         Generate_params.PriorMean = mean(log(Generate_params.ratings(:,num_subs_found)));
        %         Generate_params.PriorVar = var(log(Generate_params.ratings(:,num_subs_found)));
        %             Generate_params.BVrange = log(Generate_params.rating_bounds);   %Used for normalising BV
        %             list.allVals = log(squeeze(Generate_params.seq_vals(sequence,:,num_subs_found)));
        %         else
        
        %             Generate_params.BVrange = Generate_params.rating_bounds;    %Used for normalising BV
        
        list.allVals                = squeeze(Generate_params.seq_vals(sequence,:,num_subs_found));
        Generate_params.PriorMean   = mean(Generate_params.ratings(:,num_subs_found));
        Generate_params.PriorVar    = var(Generate_params.ratings(:,num_subs_found));
        %         end;
        
        % ranks for this sequence
        dataList                    = tiedrank(squeeze(Generate_params.seq_vals(sequence,:,num_subs_found))');
        
        %         list.optimize     = 0;
        %         list.flip         = 1;
        list.vals                   =  list.allVals;
        %         list.length       = Generate_params.seq_length;
        
        
    end % end of sequences loop
    
    
end % end of subjects_found loop



return