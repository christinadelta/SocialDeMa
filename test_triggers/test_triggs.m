%% ---------------------------------------

% Test biosemi triggers using a photosensor

% Created to test the triggers of the Optimal Stopping Problems Project
% using a BIOSEMI equipment and ActiView acquisition Software 

% Depependencies:
% 1. Matlab 2019a 
% 2. Psychtoolbox 3

% CLEAN UP
clear;
clc
close all hidden;

%% ---------------------------------------
% DEFINE INIT VARIABLES 

fixation    = '+';      % fixation cross 
EEG         = 0;        % set to 1 when running in the eeglab
triggerdur  = 0.0003;   % duration of the trigger in sec (3 ms) 
jitter      = .4;       % 0.4 sec
isi         = .5;       % in seconds
stim_dur    = 1;        % bead duration in seconds
fix_dur     = 0.5;      % duration of the fixation cross
visangle    = 4;        % degrees of visual angle

trigger1    = 1;
trigger2    = 2;

%% ---------------------------------------
% INITIAL SETUP 

basedir         = pwd;
stimdir         = fullfile(basedir, 'stimuli');

addpath(genpath(fullfile('stimdir'))); 

% Add PTB to your path and start the experiment 
ptbdir          = '/Applications/Psychtoolbox';                             
addpath(genpath(ptbdir))

% how many lines of no interest do we have in the excel file? (useful to
% remove headers)
headers                 = 1; 
% how many columns?
columns                 = 2; 

% read the excel file
[vars, txt ,~]          = xlsread('test_stimuli.xls');

vars                    = vars';
% remove the headers
txt(headers,:)          = [];

data                    = [];
objects                 = length(vars);

% load the stimuli
for i=1:objects

    img             = fullfile(stimdir,txt{i});
    image           = imread(img);
    data(i).file    = image; % should resize or not?
end

try
    %% ---------------------------------------
    % PREP EXPERIMENT(open screen, etc..) 
    
    % define colours
    black              = [0 0 0];
    white              = [255 255 255];
    grey               = [128 128 128];
    
    % text settings
    textfont           = 'Verdana';
    textsize           = 22;
    fixationsize       = 30;
    
    % Screen('Preference', 'SkipSyncTests', 0) % set a Psychtoolbox global preference.
    Screen('Preference', 'SkipSyncTests', 1) % for testing I have set this to 1. When running the actuall task uncomment the above

    screenNumber            = max(Screen('Screens'));
    
    [window, windrect]      = Screen('OpenWindow',screenNumber,grey); % open window
    
    AssertOpenGL;                                                           % Break and issue an error message if PTB is not based on OpenGL or Screen() is not working properly.
    Screen('Preference', 'Enable3DGraphics', 1);                            % enable 3d graphics
    Screen('BlendFunction', window, GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);   % Turn on blendfunction for for the screen
    priorityLevel = MaxPriority(window);                                    % query the maximum priority level
    Priority(priorityLevel);
    HideCursor;
    
    [xcenter, ycenter]      = RectCenter(windrect);                         % get the centre coordinate of the window in pixels
    [xpixels, ypixels]      = Screen('WindowSize', window);                 % size of the on-screen window in pixels
    globalrect              = Screen('Rect', screenNumber);                 % this is used for the slider

    
    % pc actual screen settings
    actscreen               = Screen('Resolution', screenNumber);
    [actwidth, actheight]   = Screen('DisplaySize', screenNumber);
    acthz                   = Screen('FrameRate', window, screenNumber);    % maximum speed at which you can flip the screen buffers, we normally use the flip interval (ifi), but better store it 
    
    ifi                     = Screen('GetFlipInterval', window);            % frame duration, inverse of frame rate, returns the duration of the frame in miliseconds
    
    slack                   = Screen('GetFlipInterval', window)/2;          % Returns an estimate of the monitor flip interval for the specified onscreen window (this is frame duration /2)
    
    %% ---------------------------------------
    % CREATE TRIAL STRUCTURE
    
    reps            = 14;
    trialArray      = repmat(vars,1,reps);
    array           = length(trialArray);
    
    % CREATE FIXATION WINDOW
    fixationdisplay = Screen('OpenOffscreenWindow',window);
    Screen('FillRect', fixationdisplay, grey);
    Screen('TextFont',fixationdisplay, textfont);
    Screen('TextSize',fixationdisplay, fixationsize);
    DrawFormattedText(fixationdisplay, fixation, 'center', ycenter, white);
    
    
    %% ---------------------------------------
    % ADD THE TRIGGER INFORMATION (IF EEG = 1) 

    if EEG == 1

        % INIT COMMUNICATION WITH EXTERNAL DEVICES
        ioObj           = io64;             % create an instance of the io64 object
        status          = io64(ioObj);      % initialize the interface to the inpoutx64 system driver
        address         = hex2dec('E010');  % LPT3 output port address for windows 10 os
        
        fprintf(' >>> OPENING USB TRIGGER LINK  <<<')
    end
    
    %% ---------------------------------------
    % MAKE IMAGE TEXTURES 
    
    % make a cell to store the image textures
    textures        = cell(1,objects);
    
    for i=1:objects
        
        % make textures
        textures{i} = Screen('MakeTexture', window, data(i).file); 
        
    end
    
    %% ---------------------------------------
    % BEGIN TESTING
   
    % prepare a general window 
    generalwindow = Screen('OpenOffscreenWindow', window, windrect);
    Screen('TextSize', generalwindow, textsize);
    Screen('FillRect', generalwindow, grey ,windrect);
    DrawFormattedText(window,'Testing will start shortly','center',ycenter,white);
    expstart = Screen('Flip', window);
    object_offset = expstart + 3;
    
    % START THE TESTING WITH A FIXATION CROSS
    Screen('CopyWindow', fixationdisplay,window, windrect, windrect)
    object_onset    = Screen('Flip', window, object_offset - slack); % flip fixation window
    
    % object offset
    object_offset   = object_onset + isi + randperm(jitter*1000,1)/1000 - ifi;
    
    for trial = 1:array
        
        thisitem    = trialArray(1,trial);
        stimtrigger = thisitem;
        
        % display the image at requested onset
%         Screen('DrawTexture', window, textures{thisitem}, [], destrect);
        Screen('DrawTexture', window, textures{thisitem});
        object_onset        = Screen('Flip', window, object_offset - slack);    % flip window
        
        % send stimulus trigger
        if EEG == 1 
            io64(ioObj, address, stimtrigger)
            WaitSecs(triggerdur);
            io64(ioObj, address, 0) % return port to zero
        end
        
        object_offset       = object_onset + stim_dur - ifi;
          
    end % end of trial loop
    
    WaitSecs(2); % wait two sec 
    
    %% ---------------------------------------
    % END OF TESTING
    
    % THIS IS IT...
    % show thank you window
    Screen('OpenOffscreenWindow', window, windrect);
    Screen('TextSize', window, textsize);
    Screen('TextStyle', window, 0)
    Screen('FillRect', window, grey ,windrect);
    DrawFormattedText(window, 'This is the end of the triggers testing.', 'center', ycenter, white);
    Screen('Flip',window, object_offset - slack);
    WaitSecs(3);
    
    % clean up at the end of testing
    Screen('CloseAll');
    ShowCursor;
    Priority(0);
    fclose('all');
   
catch % catch last errors
    
    Screen('CloseAll');
    ShowCursor;
    Priority(0);
    psychrethrow(psychlasterror);
    
end
