% PRE-PROCESSING SCRIPT FOR BEADS TASK

% Part of the Optimal Stopping Problems Project

%% IMPORTANT NOTE %%

% I save two different types of mat files. In 4 mat files I save the block
% information. This file contains the the trials of the given block, the
% number of draws for each trial, the balance of this trial (reward/loss),
% response, accuracy, etc..

% The rest mat files contain sequence/trial information, such as: each draw
% of the given trial, trial start, bead onset, response time, etc...

% TOTAL MAT FILES for each subject: 56
% BLOCK MAT FILES: 4
% SEQUENCE MAT FILES: 52

% I first extract block data, remove nans and save the data in a csv file 
% Then I extract the sequence data, remove nans and save the data in a csv file 

%% LOAD DATA %%

% GET PATHS & DEFINE PATHS
startpath       = '/Users/christinadelta/githubstuff/rhul_stuff/SocialDeMa/';
resultspath     = fullfile(startpath, 'experiments', 'results');
task            = 'beads';
subpath         = fullfile(startpath, task);
session         = 1;

subs            = dir(fullfile(resultspath, task, '*sub*'));
nsubs           = length(subs);

totaltrials     = 52; 
blocktrials     = 13;
blocks          = 4;
conditions      = 2;
thisequence     = nan(1,blocktrials); % every sequence has different number of draws
temp            = 0;

% only keep subnames
subname         = {subs.name};

%% EXTRACT AND SAVE THE BLOCK DATA %%

for subI = 1:nsubs 
    
    fprintf('loading beads block data\n')  
    subject = subs(subI).name;
    subdir = fullfile(resultspath, task,subject);
    fprintf('\t reading data from subject %d\n',subI); 
    
    for blockI = 1:blocks
        
        fprintf('\t\t loading block %d\n\n',blockI);
        subFile = fullfile(subdir, sprintf('subject_%02d_task_%s_block_%02d_ses_%02d_logs.mat',subI, task, blockI,session));
        load(subFile)
        
        for trial = 1:blocktrials
            
            indx = ((blockI -1)*blocktrials) + trial;  
            
            block(indx)     = logs.blocktrials(trial).block;
            trialno(indx)   = logs.blocktrials(trial).trialnumber;
            urntype(indx)   = logs.blocktrials(trial).urntype;
            draws(indx)     = logs.blocktrials(trial).draws;
            response(indx)  = logs.blocktrials(trial).response;
            accuracy(indx)  = logs.blocktrials(trial).accuracy;
            rate(indx)      = logs.blocktrials(trial).thisrate;
            condition(indx) = logs.blocktrials(trial).condition;
            balance(indx)   = logs.blocktrials(trial).balance;
            sub(indx)       = subI;
            
        end % end of trial loop
        
    end % end of block loop
    
end % end of subject loop

% add data in one matrix
data = [sub' block' trialno' urntype' draws' response' accuracy' rate' condition' balance'];

% remove nans if any
data(any(isnan(data), 2), :)  = [];

% % save matrix in csv format for r and python
csvwrite('beads_data.csv', data)




% % load subject
% for sub = 1:nsubs
%     
%     fprintf('loading beads block data\n')  
%     subject = subs(sub).name;
%     subdir = fullfile(resultspath, task,subject);
%     fprintf('\t reading data from subject %d\n',sub); 
%     
%     % load block
%     for blockI = 1:blocks
%         
%         fprintf('\t\t loading block %d\n\n',block);
%         
%         % load trial
%         for trial = 1:blocktrials
%             
%             % load sequence file
%             fprintf('\t\t loading sequence %d\n\n',trial);
%             subFile = fullfile(subdir, sprintf('subject_%02d_task_%s_block_%02d_trial_%02d_ses_%02d_draw_logs.mat',sub, task, blockI, trial,session));
%             load(subFile)                     
%             
%             % how many draws in that sequence?
%             thisequence(trial) = length(logs.draws);
%             
%             for i = 1:thisequence(trial)
%                 
%                 index               = temp + i;    
% 
%                 block(index)        = logs.draws(i).block;
%                 trialno(index)      = logs.draws(i).trialnumber;
%                 drawno(index)       = logs.draws(i).thisdraw;
%                 bead(index)         = logs.draws(i).thisbead;
%                 rt(index)           = logs.draws(i).rt;
%                 
%             end
%             temp = temp + thisequence(trial); % update temp 
%         end
%     end
%     
% end % end of subject loop
