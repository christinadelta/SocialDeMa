%% ---------------------------------------

%%%%%%%%%%%%%% BEADS TASK V3 %%%%%%%%%%%%%

% DESCRIPTION:

% Third version of the Beads task implemented with the PsychToolbox
% Depependencies:
% 1. Matlab 2021a 
% 2. Psychtoolbox 3

% For helpful info regarding the psychtoolbox see:
% http://peterscarfe.com/ptbtutorials.html

%%% TODO: %%%
% add logs txt file? 
% now we save only in .mat files 

% CLEAN UP
% clear;
% clc
% close all hidden;

%% ---------------------------------------
% INITIAL EXPERIMENTAL SETUP 

% Initialize the random number generator
rand('state', sum(100*clock)); 


% get participant nb and task name 
answer          = startup.answer;

% initial experimental settings
sub             = str2num(answer{2}); % participant number
taskName        = answer{1}; 
taskNb          = 1; 
sess            = 1;

basedir         = pwd;

% get directories and add utility functions to the path
wd      = fullfile(basedir, 'experiments');
addpath(genpath(fullfile(wd,'utils')));                              % add subfunctions to the path

%% ---------------------------------------
% SET OUTPUT INFO AND LOGS FILE

logs.sub            = sub;
logs.task           = taskName;
logs.sess           = sess;
logs.date           = datestr(now, 'ddmmyy');
logs.time           = datestr(now, 'hhmm');

logs.drawslog       = 'subject_%02d_task_%s_block_%02d_trial_%02d_ses_%02d_draw_logs.mat';
logs.trialog        = 'subject_%02d_task_%s_block_%02d_ses_%02d_logs.mat';
logs.txtlog         = 'subject_%02d_task_%s_block_%02d_trial_%02d_ses_%02d_events.tsv';

% % setup study output file
logs.resultsfolder  = fullfile(wd, 'results',taskName, sprintf('sub-%02d', sub));

if ~exist(logs.resultsfolder, 'dir')
    mkdir(logs.resultsfolder)
end

% Add PTB to your path and start the experiment 
ptbdir          = '/Applications/Psychtoolbox';                             % change to your ptb directory
addpath(genpath(ptbdir))

scrn.ptbdir     = ptbdir;

try
    %% ---------------------------------------
    % PREP EXPERIMENT(open screen, etc..) 
    
    % define colours
    scrn.black              = [0 0 0];
    scrn.white              = [255 255 255];
    scrn.grey               = [128 128 128];
    scrn.green              = [0 140 54];
    scrn.blue               = [30 70 155];
    scrn.red                = [225 25 0];
    scrn.orange             = [243 146 0];

    % text settings
    scrn.textfont           = 'Verdana';
    scrn.textsize           = 20;
    scrn.smalltext          = 12;
    scrn.fixationsize       = 30;
    scrn.textbold           = 1; 
    
    % Screen('Preference', 'SkipSyncTests', 0) % set a Psychtoolbox global preference.
    Screen('Preference', 'SkipSyncTests', 1) % for testing I have set this to 1. When running the actuall task uncomment the above

    screenNumber            = max(Screen('Screens'));
    
    [window, windrect]      = Screen('OpenWindow',screenNumber, scrn.grey); % open window
    
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
    scrn.acthz              = Screen('FrameRate', window, screenNumber);    % maximum speed at which you can flip the screen buffers, we normally use the flip interval (ifi), but better store it 
    
    scrn.ifi                = Screen('GetFlipInterval', window);            % frame duration, inverse of frame rate, returns the duration of the frame in miliseconds
    
    scrn.slack              = Screen('GetFlipInterval', window)/2;          % Returns an estimate of the monitor flip interval for the specified onscreen window (this is frame duration /2)
    
    % create rectangle for trials
    screenresolution        = [actscreen.width actscreen.height];
    baserect                = [0 0 ceil(.5*screenresolution(1)) ceil(0.5*screenresolution(2))];
    %baserect                = [0 0 800 800];
    centeredrect            = CenterRectOnPointd(baserect, xcenter, ycenter);
    penwidth                = 4; % line width
    
    % UPDATE SCREEN STRUCT
    scrn.actscreen          = actscreen;
    scrn.frame_rate         = 1/scrn.ifi;
    scrn.actwidth           = actwidth;
    scrn.actheight          = actheight;
    scrn.window             = window;
    scrn.windrect           = windrect;
    scrn.xcenter            = xcenter;
    scrn.ycenter            = ycenter;
    scrn.xpixels            = xpixels;
    scrn.ypixels            = ypixels;
    scrn.globalrect         = globalrect;
    scrn.screenNumber       = screenNumber;
    scrn.screenresolution   = screenresolution;
    scrn.baserect           = baserect;
    scrn.centeredrect       = centeredrect;
    scrn.penwidth           = penwidth;
    
    %% ---------------------------------------
    % RUN A FEW IMPORTANT UTIL FUNCTIONS
    
    set                 = TaskSettings(taskNb, sess);                       % Define the first task-specific parameters

    set                 = DefineKeys(taskNb, set);                          % Define set of the task
    
    scrn                = screenSettings(scrn, taskNb);                     % Define screen setup
    
    [trials, set]       = CreateTrialList(set);                             % create trials, sequences, split in runs, etc..
    
    %% ---------------------------------------
    % CREATE AND RUN INSTRUCTIONS
    
    % UNPACK SETTINGS
    iduration           = set.welcomedur;
    EEG                 = set.EEG; % should be EEG = 1 when running at the EEGlab
    spacekey            = set.code20;
    esckey              = set.code21;
    
    % prepare a general window (this will be used for instruction and
    % information display 
    generalwindow = Screen('OpenOffscreenWindow', window, windrect);
    Screen('TextSize', generalwindow, scrn.textsize);
    Screen('FillRect', generalwindow, scrn.grey ,windrect);
    
    % Start instructions
    DrawFormattedText(window,'Hello! Please pay attention to the instructions','center',scrn.ycenter,scrn.white);
    expstart = Screen('Flip', window);
    duration = expstart + iduration;
   
    
    % display instructions window 1
    instructions = Screen('OpenOffscreenWindow', window, windrect);
    Screen('TextSize', instructions, scrn.textsize);
    Screen('FillRect', instructions, scrn.grey ,windrect);
    DrawFormattedText(instructions, 'There are two urns:', 'center', scrn.ycenter-250, scrn.white);
    DrawFormattedText(instructions, 'The Blue Urn has more blue balls than green balls. The Green Urn has more green balls than blue balls.', 'center', scrn.ycenter-200, scrn.white);
    DrawFormattedText(instructions, 'On each trial, you will draw a sequence of balls from one of these two urns. Your job is to decide whether', 'center', scrn.ycenter-150, scrn.white);
    DrawFormattedText(instructions, 'the balls are drawn from the blue urn or the green urn. After each ball is drawn, you may choose to: ','center', scrn.ycenter-100, scrn.white);
    DrawFormattedText(instructions, 'Guess The Blue Urn by pressing the keycode 1', 'center', scrn.ycenter-50, scrn.white);
    DrawFormattedText(instructions, 'Guess The Green Urn by pressing the keycode 2', 'center', scrn.ycenter, scrn.white); 
    DrawFormattedText(instructions, 'Draw another ball by pressing the keycode 3', 'center', scrn.ycenter+50, scrn.white);
    DrawFormattedText(instructions, 'You may make a decision after any draw but you may not draw more than 9 balls', 'center', scrn.ycenter+100, scrn.white);
    DrawFormattedText(instructions, 'If you have understood the instructions so far, press SPACE to change page', 'center', scrn.ycenter+150, scrn.white); 
    
    % copy the instructions window  and flip.
    Screen('CopyWindow',instructions,window,windrect, windrect);
    Screen('Flip', window, duration);
    
    % WAIT FOR THEM TO PRESS SPACE
    responsemade = 1;
    while responsemade
        [~, secs, keycode]  = KbCheck;
        WaitSecs(0.001) % delay to prevent CPU logging

        % spacebar is pressed 
        if keycode(1, spacekey)
            responsemade    = 0;
        end
    end
    
    WaitSecs(1); % wait one sec before flipping to the block/trial/sequence information screen 
    
    % display instructions window 2
    instructions2 = Screen('OpenOffscreenWindow', window, windrect);
    Screen('TextSize', instructions2, scrn.textsize);
    Screen('FillRect', instructions2, scrn.grey ,windrect);
    DrawFormattedText(instructions2, 'After an urn is chosen, you will be asked to rate how confident you are about', 'center', scrn.ycenter-200, scrn.white);
    DrawFormattedText(instructions2, 'the choice that you made on a scale of 0 to 100.', 'center', scrn.ycenter-150, scrn.white);
    DrawFormattedText(instructions2, 'You will be asked to give your rating using a slider scale. When the scale appears, ', 'center', scrn.ycenter-100, scrn.white);
    DrawFormattedText(instructions2, 'you will first need to (left) click the mouse one time. This will reveal','center', scrn.ycenter-50, scrn.white);
    DrawFormattedText(instructions2, 'a black vertical line (the slider). Scroll the slider using the mouse (from left to right) and stop', 'center', scrn.ycenter, scrn.white);
    DrawFormattedText(instructions2, 'and click on the desired rating/position. Before starting the actual experiment,', 'center', scrn.ycenter+50, scrn.white)
    DrawFormattedText(instructions2, 'you will first complete a short and simple quiz, to test your understanding of the task instructions.', 'center', scrn.ycenter+100, scrn.white)
    DrawFormattedText(instructions2, 'If you have understood the instructions, press SPACE to start the short quiz.', 'center', scrn.ycenter+150, scrn.white)
    
    % copy the instructions2 window  and flip.
    Screen('CopyWindow',instructions2,window,windrect, windrect);
    Screen('Flip', window, duration);
    
    % WAIT FOR THEM TO PRESS SPACE
    responsemade = 1;
    while responsemade
        [~, secs, keycode]= KbCheck;
        WaitSecs(0.001) % delay to prevent CPU logging

        % spacebar is pressed 
        if keycode(1, spacekey)
            responsemade = 0;
        end
    end
    
    WaitSecs(1); % wait one sec before flipping to the Instructions quiz 
    
    %% ---------------------------------------
    % RUN THE INSTRUCTIONS QUIZ 
    
%     % Start instructions
%     DrawFormattedText(window,'INSTRUCTIONS QUIZ','center',scrn.ycenter,scrn.white);
%     Screen('Flip', window);
%     WaitSecs(1);
%     
%     set = ShortQuiz(set, scrn, set); % RUN the instructions quiz 
%     

    %% ---------------------------------------
    % ADD THE TRIGGER INFORMATION (IF EEG = 1) 

    if EEG == 1

        % INIT COMMUNICATION WITH EXTERNAL DEVICES
        sp      = BioSemiSerialPort();
        set.sp  = sp;
        fprintf(' >>> OPENING USB TRIGGER LINK  <<<')

        % UNPACK TRIGGERS FROM THE SETTINGS STRUCTURE
        trigger100 = set.trigger100; % condition trigger -- condition (easy)
        trigger101 = set.trigger101; % condition trigger -- condition (difficult) 
        trigger102 = set.trigger102; % start of sequence 
        trigger103 = set.trigger103; % end of sequence
    end

    %% ---------------------------------------
    % START THE BLOCK & SEQUENCE/TRIAL LOOPS
    
    abort           = 0;                % when 1 subject can quit the experiment
    
    % UNPACK SETTINGS STRUCT
    ntrials         = set.trials;       % total trials
    nb_blocks       = set.blocks;       % total blocks
    trialsPerBlock  = set.blocktrials;  % trials per block
    
    % INIT BLOCKS LOOP
    for iBlock = 1:nb_blocks
        
        % init block-log struct
        blocktrials     = [];
        
        % first allow subject to exit experiment if they pressed the esc key 
        [keyisdown,secs,keycode] = KbCheck;
        if keyisdown && keycode(esckey)
            
            % if the subject pressed ESC
            responsemade = 1;
            while responsemade
                [~, secs, keycode] = KbCheck;
                WaitSecs(0.001) % delay to prevent CPU logging

                % ESC is pressed 
                if keycode(1, esckey)
                    abort           = 1;
                    responsemade    = 0;
                end
            end
        end
        if abort == 1
            break;
        end
        
        % UNPACK TRIALS STRUCT
        block_seq           = trials.sequence{iBlock};
        block_urns          = trials.urns{iBlock};
        
        % INIT TRIALS LOOP
        for thistrial = 1:trialsPerBlock
            
            % show the sequence information screen and wait until subject
            % presses space to continue 
            
            set.thisblock   = iBlock;                       % send the current block number to the RUN fuction
            set.thistrial   = thistrial;                    % send the current trial number to the RUN function
            set.sequence    = block_seq{thistrial};         % send the current sequence to the RUN function
            set.urn         = block_urns(thistrial);        % send the current urn to the RUN function
            
            currentbalance  = set.balance;                  % show current balance at the beggining og each trial
            loss            = set.loss;
            % count the proportions of bead colours in the current sequence
            % list to dettermine the difficulty condition
            condition       = sum(set.sequence(1,:)==2);
            
            if condition == 2 % if 2 appears 2 times in the sequence
                
                cond        = 1; % add this to the blocktrials struct
                 
                high_p      = 80;
                low_p       = 20;
            else % if 2 appears 4 times in the sequence 
                
                cond        = 2;
                high_p      = 60;
                low_p       = 40;
            end
   
            % display trial/sequence information window 
            Screen('OpenOffscreenWindow', window, windrect);
            Screen('TextSize', window, scrn.textsize);
            Screen('TextStyle', window, 0)
            Screen('FillRect', window, scrn.grey ,windrect);
            DrawFormattedText(window, sprintf('Starting sequence %d of block %d',thistrial, iBlock), 'center', scrn.ycenter-100, scrn.white);
            DrawFormattedText(window, sprintf('Your current balance is: %3.3f credits',currentbalance), 'center', scrn.ycenter-50, scrn.white);
            DrawFormattedText(window, sprintf('The urns have a %02d:%02d color split. You will lose %d credits if you are wrong',high_p, low_p, loss), 'center', scrn.ycenter, scrn.white);
            DrawFormattedText(window, 'Press SPACE to continue, or press ESC to quit', 'center', scrn.ycenter+50, scrn.white);
            Screen('Flip', window); 
            
            % send a condition trigger
            if EEG == 1
                if condition == 2
                    sp.sendTrigger(trigger100)
                else
                    sp.sendTrigger(trigger101)
                end
            end
            
            % WAIT FOR THEM TO PRESS SPACE
            responsemade = 1;
            while responsemade
                [~, secs, keycode]= KbCheck;
                WaitSecs(0.001) % delay to prevent CPU logging

                % spacebar is pressed 
                if keycode(1, spacekey)
                    responsemade    = 0;
                    
                    % or esc is pressed
                elseif keycode(1, esckey)
                    abort           = 1;
                    responsemade    = 0;
                end
            end
            if abort == 1 % exit 
                break;
            end
            
            [set,logs]     = RunBeads(set, scrn, logs);
            
            % UNPACK SET AND ADD THE TRIAL INFO TO THE "BLOCK"-LOG FILE 
            blocktrials(thistrial).session     = set.blocktrials.session;
            blocktrials(thistrial).block       = set.blocktrials.block;
            blocktrials(thistrial).trialnumber = set.blocktrials.trialnumber;
            blocktrials(thistrial).trialonset  = set.blocktrials.trialonset;
            blocktrials(thistrial).urntype     = set.blocktrials.urntype;
            blocktrials(thistrial).sequence    = set.blocktrials.sequence;
            blocktrials(thistrial).draws       = set.blocktrials.draws;
            blocktrials(thistrial).response    = set.blocktrials.response;
            blocktrials(thistrial).accuracy    = set.blocktrials.accuracy;
            blocktrials(thistrial).balance     = set.blocktrials.balance;
            blocktrials(thistrial).thisrate    = set.blocktrials.thisrate;
            blocktrials(thistrial).ratert      = set.blocktrials.ratert;
            blocktrials(thistrial).condition   = cond;
               
        end % end of trial/sequence for loop        
        
        % save trial info
        logs.blocktrials    = blocktrials;
        sub_log             = fullfile(logs.resultsfolder,sprintf(logs.trialog,sub,taskName,iBlock,sess));
        save(sub_log,'logs');

        % IF THIS IS THE LAST BLOCK BREAK FROM THE LOOP AND GO DIRECTLY TO
        % THE GOODBYE SCREEN
        if iBlock == nb_blocks
            break;
        end
        
        % if this is the end of the block 
        % allow subject to take a short break (if they want to)
        Screen('OpenOffscreenWindow', window, windrect);
        Screen('TextSize', window, scrn.textsize);
        Screen('TextStyle', window, 0)
        Screen('FillRect', window, scrn.grey ,windrect);
        DrawFormattedText(window, 'Time for a break! When ready to continue, press SPACE', 'center', scrn.ycenter-50, scrn.white);
        DrawFormattedText(window, 'Or press ESC to quit.', 'center', scrn.ycenter, scrn.white);
        Screen('Flip', window); 
        
        % WAIT FOR THEM TO PRESS SPACE
        responsemade = 1;
        while responsemade
            [~, secs, keycode]= KbCheck;
            WaitSecs(0.001) % delay to prevent CPU logging

            % spacebar is pressed 
            if keycode(1, spacekey)
                responsemade    = 0;
                
            % or esc is pressed
            elseif keycode(1, esckey)
                abort           = 1;
                responsemade    = 0;
            end
        end
        
        WaitSecs(1) % wait one sec before moving to the next window
        
        if abort == 1 % exit if subject pressed ESC
            break;
        end
        
    end % end of block for loop
    
    currentbalance      = set.balance;                  % this is your balance in credits  
    conversion          = set.conversion;               % conversion rate
    totals              = currentbalance * conversion;  % this is your converted winnings
    totalreward         = totals * 0.03;                % actual winnings 
    set.totalreward     = totalreward;
        
    % THIS IS IT...
    % show thank you window
    Screen('OpenOffscreenWindow', window, windrect);
    Screen('TextSize', window, scrn.textsize);
    Screen('TextStyle', window, 0)
    Screen('FillRect', window, scrn.grey ,windrect);
    DrawFormattedText(window, 'This is the end of the experiment.', 'center', scrn.ycenter-50, scrn.white);
    DrawFormattedText(window, sprintf('Your calculated total reward is: £%3.3f', totalreward), 'center', scrn.ycenter, scrn.white);
    DrawFormattedText(window, 'Thank you for your time!', 'center', scrn.ycenter+50, scrn.white);
    Screen('Flip',window);
    WaitSecs(3);
    
    % clean up at the end of the experiment
    Screen('CloseAll');
    ShowCursor;
    Priority(0);
    fclose('all');

catch % catch last errors
    
    Screen('CloseAll');
    ShowCursor;
    Priority(0);
    psychrethrow(psychlasterror);
    
end % end try... catch