function [set] = loaditems(set, wd)

% SUBFUNCTION USED IN:
% Best Choice Economic Task
% Best Choice Facial Attractiveness Task

% this subfucntion load the items from the excel files of the best-coice
% tasks 

taskNb              = set.taskNb;
exceldir            = fullfile(wd, 'excel_files');

if taskNb == 2
    
    % how many lines of no interest do we have in the excel file? (useful to
    % remove headers)
    headers                 = 1; 
    % how many columns?
    columns                 = 7; 
    
    % read the excel file
    [vars, txt ,~]          = xlsread(fullfile(exceldir, 'BestChoice_Economic_ItemList.xlsx'));
    
    % remove the headers
    txt(headers,:)          = [];
    
    set.items               = vars;

    % keep the relevant information
    for i = 1:length(vars)
        set.data{i}        = txt{i,5};
        set.description{i} = txt{i,4};
        set.price{i}       = txt{i,7};
        set.model{i}       = txt{i,8};
    end
    
end % end of taskNb if statement 


return