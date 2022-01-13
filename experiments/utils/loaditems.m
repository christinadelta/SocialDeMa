function [set] = loaditems(set, wd)

% SUBFUNCTION USED IN:
% Best Choice Economic Task
% Best Choice Facial Attractiveness Task

% this subfucntion load the items from the excel files of the best-coice
% tasks 

taskNb              = set.taskNb;
exceldir            = fullfile(wd, 'excel_files');

% which task is it? which stimuli to load?
if taskNb == 3
    type = set.type;
end

if taskNb == 2
    
    % how many lines of no interest do we have in the excel file? (useful to
    % remove headers)
    headers                 = 1; 
    % how many columns?
    columns                 = 5; 
    
    % read the excel file
    [vars, txt ,~]          = xlsread(fullfile(exceldir, 'economic_best_choice.xls'));
    
    % remove the headers
    txt(headers,:)          = [];
    vars(headers,:)         = [];
    
    set.items               = vars(:,1);
    set.price               = vars(:,4);

%     % keep the relevant information
%     for i = 1:length(vars)
%         set.data{i}        = txt{i,3};
%         set.description{i} = txt{i,2};
%         set.price{i}       = vars(i,4);
%         set.model{i}       = txt{i,5};
%     end
    
elseif taskNb == 3
    
    imgdir                  = fullfile(wd, 'stimuli');
    
    % male or female stimuli?
    if type == 1
        stimdir             = fullfile(imgdir, 'female');
    else
        stimdir             = fullfile(imgdir, 'male');
    end
    
    % how many lines of no interest do we have in the excel file? (useful to
    % remove headers)
    headers                 = 1; 
    % how many columns?
    columns                 = 2; 
    
    % read the excel file
    [vars, txt ,~]          = xlsread(fullfile(exceldir, 'face_filenames.xls'));
    
    % remove the headers
    txt(headers,:)          = [];
    
    data                    = [];
    objects                 = length(vars);
     
    % load the stimuli
    for i=1:objects
        
        img             = fullfile(stimdir,txt{i});
        image           = imread(img);
        data(i).file    = imresize(image,[set.stimsize set.stimsize]); % should resize or not?
    end
    
    stimwidth           = size(data(1).file,1);   % width of objects
    stimhight           = stimwidth;
    
    % UPDATE SETTINGS STRUCT
    set.stimwidth       = stimwidth;
    set.stimhight       = stimhight;
    set.objects         = objects;
    set.data            = data;
    set.items           = vars; 
    
end % end of taskNb if statement 


return