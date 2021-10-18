% PRE-PROCESSING SCRIPT FOR BEST-CHOICE ECONOMIC TASK

% get paths and define vars
startpath       = '/Users/christinadelta/githubstuff/rhul_stuff/SocialDeMa/';
resultspath     = fullfile(startpath, 'experiments', 'results');
task            = 'economic';
subpath         = fullfile(startpath, task);
phases          = 2;
session         = 2;

subs            = dir(fullfile(resultspath, task, '*sub*'));
nsubs           = length(subs);

totaltrials     = 40; 
blocktrials     = 20;
blocks          = 2;
phase           = 2;

% only keep subnames
subname         = {subs.name};

for subI = 1:nsubs 
    
    fprintf('loading economic best-choice data\n')  
    subject = subs(subI).name;
    subdir = fullfile(resultspath, task,subject);
    fprintf('\t reading data from subject %d\n',subI);  
    
    for blockI = 1:blocks
        
        fprintf('\t\t loading block %d\n\n',blockI);
        subFile = fullfile(subdir, sprintf('subject_%02d_task_%s_block_%02d_ses_%02d_phase_%02d_blocktrials_logs.mat',subI, task, blockI,session,phase));
        load(subFile)
        
        for trial = 1:blocktrials
            
            indx = ((blockI -1)*blocktrials) + trial;  % trials in total (1440)
            
            trialno(indx) = logs.blocktrials(trial).trialnumber;
            samples(indx) = logs.blocktrials(trial).numsamples;
            item(indx) = logs.blocktrials(trial).chosenitem;
            price(indx) = logs.blocktrials(trial).chosenprice;
            rank(indx) = logs.blocktrials(trial).rank;
            reward(indx) = logs.blocktrials(trial).reward;
            balance(indx) = logs.blocktrials(trial).balance;
            sub(indx) = subI;
        end
  
    end % end of block loop
   
end % end of subjects loop

% add data in one matrix
data = [trialno' item' samples' price' rank' reward' balance'];

% % save matrix in csv format for r and python
csvwrite('economic_phase2_data.csv', data)
