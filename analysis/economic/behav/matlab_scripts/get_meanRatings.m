function [phase1_data, meanRatings] = get_meanRatings(subdir,sub,task)

% This function is used with Phase 1 -- Economic data

session     = 1;
phase       = 1;
blocktrials = 40;
blocks      = 20;
allitems    = 400;

% we have to subjects with incomplete datasets :/
if sub == 7 
    blocks = 17;
elseif sub == 40 
    blocks = 15;
end

%% EXTRACT DATA %%

for block = 1:blocks

    fprintf('\t\t loading block %d\n\n',block);
    subFile = fullfile(subdir, sprintf('subject_%02d_task_%s_block_%02d_ses_%02d_phase_%02d_logs.mat',sub, task, block,session,phase));
    load(subFile)

    % extract trial info
    for trial = 1:blocktrials

        indx            = ((block -1)*blocktrials) + trial;  

        subj(indx)      = sub;
        trialno(indx)   = logs.trials(trial).trialNb;
        blockno(indx)   = logs.trials(trial).block;
        thisitem(indx)  = logs.trials(trial).thisitem;
        thisprice(indx) = logs.trials(trial).thisprice;
        rate(indx)      = logs.trials(trial).response;
        rt(indx)        = logs.trials(trial).rt;


    end % end of trials loop


end % end of blocks loop

% store all in one matrix
phase1_data = [subj' blockno' trialno' thisitem' thisprice' rate' rt'];

clear subj trialno blockno thisitem thisprice rate rt session phase indx
%%  GET AVERAGE OF EACH RATING %% 

% first, how many unique prices where there to rate?
uitems                      = length(unique(phase1_data(:,4)));


% init averaged ratings array
ratings                 = nan(uitems,1);
items                   = [1:uitems]'; % array [1:400]

% loop over unique prices and average rating for each price
for i = 1:allitems 
    
    tmpitem             = find(phase1_data(:,4) == i); % find this_item

    if isempty(tmpitem) % this is for incomplete datasets
        tempitem = i; % 
        tmprate             = nan;
        tmpprice            = nan;
    else
        tmprate             = phase1_data((tmpitem),6);
        tmpprice            = phase1_data((tmpitem),5);
    end

   
    % average prices for item i and store
    ratings(i,1)        = i;                % store item number
    ratings(i,2)        = nanmean(tmprate);    % store rating for that item
    ratings(i,3)        = tmpprice(1);      % store actual price for that item

    % each price is shown twise, so the tmpprice variable contains
    % the same price twise, just pick the first one
    items(i,2)          = tmpprice(1);
    
end % end of unique prices loop

meanRatings             = ratings;


end % end of function