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
        
        list.allVals                = squeeze(Generate_params.seq_vals(sequence,:,num_subs_found)); % squeezing may not be needed
        Generate_params.PriorMean   = mean(Generate_params.ratings(:,num_subs_found));
        Generate_params.PriorVar    = var(Generate_params.ratings(:,num_subs_found));
        %         end;
        
        % ranks for this sequence
        dataList                    = tiedrank(squeeze(Generate_params.seq_vals(sequence,:,num_subs_found))');
        
        %         list.optimize     = 0;
        %         list.flip         = 1;
        list.vals                   =  list.allVals;
        %         list.length       = Generate_params.seq_length;
        
        % Do cutoff model, if needed
        if Generate_params.model(Generate_params.current_model).identifier == 1
            
            % get seq vals to process
            this_seq_vals   = list.allVals;
            % initialise all sequence positions to zero/continue (value of stopping zero)
            choiceStop      = zeros(1,Generate_params.seq_length);
            % What's the cutoff?
            estimated_cutoff = round(Generate_params.model(Generate_params.current_model).cutoff);
            if estimated_cutoff < 1; estimated_cutoff = 1; end
            if estimated_cutoff > Generate_params.seq_length; estimated_cutoff = Generate_params.seq_length; end
            
            % find seq vals greater than the max in the period
            % before cutoff and give these candidates a maximal stopping value of 1
            choiceStop(1,find( this_seq_vals > max(this_seq_vals(1:estimated_cutoff)) ) ) = 1;
            % set the last position to 1, whether it's greater than
            % the best in the learning period or not
            choiceStop(1,Generate_params.seq_length) = 1;
            % find first index that is a candidate ....
            num_samples(sequence,this_sub) = find(choiceStop == 1,1,'first');   %assign output num samples for cut off model
            
            % Reverse 0s and 1's for ChoiceCont
            choiceCont = double(~choiceStop);
            
        else % if any of the Bayesian models
            
            % run Nick's model
            [choiceStop, choiceCont, difVal]  = ...
                analyzeSecretaryNick_2021(Generate_params,list);
            
            num_samples(sequence,this_sub) = find(difVal<0,1,'first');  % assign output num samples for Bruno model
            
            
        end % end of if statement
        
        % ...and its rank
        ranks(sequence,this_sub) = dataList( num_samples(sequence,this_sub) );
        % Accumulate action values too so you can compute ll outside this function if needed
        choiceStop_all(sequence, :, this_sub) = choiceStop;
        choiceCont_all(sequence, :, this_sub) = choiceCont;
        
    end % end of sequences loop
    
    this_sub = this_sub + 1;
    
end % end of subjects_found loop



return