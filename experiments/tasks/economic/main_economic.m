%% ---------------------------------------
% DESCRIPTION:

% First version of the Economic Best-Choice task implemented with the PsychToolbox
% Depependencies:
% 1. Matlab 2021a 
% 2. Psychtoolbox 3

% For helpful info regarding the psychtoolbox see:
% http://peterscarfe.com/ptbtutorials.html

%%% TODO: %%%
% add logs txt file 
% now we save only in .mat files 

%% ---------------------------------------
% INITIAL EXPERIMENTAL SETUP 

% Initialize the random number generator
rand('state', sum(100*clock));

% get participant nb and task name 
answer          = startup.answer;

% initial experimental settings
sub             = str2num(answer{2}); % participant number
taskName        = answer{1}; 
taskNb          = 2; 

basedir         = pwd;

% add a new prompt to identify the phase of the task 
prompt          = {'Which phase is it?'};
dlgtitle        = 'Info window';
dims            = [1 30];
definput        = {'1'}; % this is a default input (this should change)
taskphase       = inputdlg(prompt,dlgtitle,dims,definput);

phase           = str2num(taskphase{1});                                    % convert task phase to number
sess            = phase;
% get directories and add utility functions to the path
wd              = fullfile(basedir, 'experiments');
addpath(genpath(fullfile(wd,'utils')));                                     % add subfunctions to the path


%% ---------------------------------------
% SET OUTPUT INFO AND LOGS FILE

logs.sub                = sub;
logs.task               = taskName;
logs.sess               = sess;
logs.date               = datestr(now, 'ddmmyy');
logs.time               = datestr(now, 'hhmm');

logs.trialog            = 'subject_%02d_task_%s_block_%02d_ses_%02d_phase_%02d_logs.mat';
logs.txtlog             = 'subject_%02d_task_%s_block_%02d_ses_%02d_phase_%02d_events.tsv';

if phase == 2
    logs.blocktrialog   = 'subject_%02d_task_%s_block_%02d_ses_%02d_phase_%02d_blocktrials_logs.mat';
end

% % setup study output file
logs.resultsfolder      = fullfile(wd, 'results',taskName, sprintf('sub-%02d', sub));

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
    scrn.black      = [0 0 0];
    scrn.white      = [255 255 255];
    scrn.grey       = [128 128 128];
    scrn.green      = [0 140 54];
    scrn.blue       = [30 70 155];
    scrn.red        = [225 25 0];

    % text settings
    scrn.textfont       = 'Verdana';
    scrn.textsize       = 20;
    scrn.fixationsize   = 30;
    scrn.textbold       = 1; 
    
    % create text settings for the previous sample (this should be very
    % small to appear at the bottom of the screen
    if phase == 2
        scrn.ptextsize      = 8;
        scrn.ptextbold      = 1;
    end
    
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

    
    % pc actual screen settings
    scrn.actscreen          = Screen('Resolution', screenNumber);
    [actwidth, actheight]   = Screen('DisplaySize', screenNumber);
    scrn.acthz              = Screen('FrameRate', window, screenNumber);    % maximum speed at which you can flip the screen buffers, we normally use the flip interval (ifi), but better store it 
    
    scrn.ifi                = Screen('GetFlipInterval', window);            % frame duration, inverse of frame rate, returns the duration of the frame in miliseconds
    
    scrn.slack              = Screen('GetFlipInterval', window)/2;          % Returns an estimate of the monitor flip interval for the specified onscreen window (this is frame duration /2)
    
    scrn.frame_rate         = 1/scrn.ifi;
    scrn.actwidth           = actwidth;
    scrn.actheight          = actheight;
    scrn.window             = window;
    scrn.windrect           = windrect;
    scrn.xcenter            = xcenter;
    scrn.ycenter            = ycenter;
    scrn.xpixels            = xpixels;
    scrn.ypixels            = ypixels;
    
    %% ---------------------------------------
    % RUN A FEW IMPORTANT UTIL FUNCTIONS
    
    set                 = TaskSettings(taskNb, sess);                       % Define the first task-specific parameters
    
    set                 = DefineKeys(taskNb, set);                          % Define keys of the task
    
    set                 = loaditems(set, wd);                               % read the excel file with the items (contracts)
    
    scrn                = screenSettings(scrn, set);                        % Define screen setup
    
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
    DrawFormattedText(window,'Hello! Please pay attentions to the instructions','center',scrn.ycenter,scrn.white);
    expstart = Screen('Flip', window);
    duration = expstart + iduration;
    
    if phase == 1 % show the instructions of phase 1 of the economic best choice task
        
        % display instructions for phase 1
        instructions = Screen('OpenOffscreenWindow', window, windrect);
        Screen('TextSize', instructions, scrn.textsize);
        Screen('FillRect', instructions, scrn.grey ,windrect);
        DrawFormattedText(instructions, 'This is the first phase of the experiment. You will be presented with smartphone contracts', 'center', scrn.ycenter-200, scrn.white);
        DrawFormattedText(instructions, 'one-by-one at the centre of the screen. Your task is to carefully view each contract and rate', 'center', scrn.ycenter-150, scrn.white);
        DrawFormattedText(instructions, '"how likely it would be to choose this contract in real life" on a scale of 1 to 9, where 1 means "I would not never choose this contract"', 'center', scrn.ycenter-100, scrn.white);
        DrawFormattedText(instructions, 'and 9 means "I would definitely choose this contract".','center', scrn.ycenter-50, scrn.white);
        DrawFormattedText(instructions, 'For you responses, press the keyboard keys 1 to 9 which correspond to your rating for each contract.', 'center', scrn.ycenter, scrn.white);
        DrawFormattedText(instructions, 'Please rate the smartphone contracts exactly as you would you in a real life scenario.', 'center', scrn.ycenter+50, scrn.white); 
        DrawFormattedText(instructions, 'If you have understood the instructions so far, press SPACE to continue', 'center', scrn.ycenter+100, scrn.white);
        
    else % if this is the second phase of the experiment 
        
        % display instructions for phase 2
        instructions = Screen('OpenOffscreenWindow', window, windrect);
        Screen('TextSize', instructions, scrn.textsize);
        Screen('FillRect', instructions, scrn.grey ,windrect);
        DrawFormattedText(instructions, 'This is the second phase of the experiment. This phase is plit in 30 sequences.', 'center', scrn.ycenter-250, scrn.white);
        DrawFormattedText(instructions, 'On each sequence, you will be presented with 10 smartphone contracts from the the previous phase, one-by-one.', 'center', scrn.ycenter-200, scrn.white);
        DrawFormattedText(instructions, 'Every time you are presented with a contract, you may either "chooce to accept it" or you may "reject it" and', 'center', scrn.ycenter-150, scrn.white);
        DrawFormattedText(instructions, 'view the next contract. Please note that for each sequence you can choose only one contract. If you reject','center', scrn.ycenter-100, scrn.white);
        DrawFormattedText(instructions, 'a contract, you may not go back and choose it. If, by the end of a sequence you have not chosen', 'center', scrn.ycenter-50, scrn.white);
        DrawFormattedText(instructions, 'a contract, by default, the last contract will be saved as your chosen contract.', 'center', scrn.ycenter, scrn.white); 
        DrawFormattedText(instructions, 'Note that each contract that you reject will be desplayed at the bottom of the screen', 'center', scrn.ycenter+50, scrn.white);
        DrawFormattedText(instructions, 'so that you have an idea of the contracts you rejected, and the number of contracts left in the sequence.', 'center', scrn.ycenter+100, scrn.white);
        DrawFormattedText(instructions, 'Press keyboard key "1" to reject a contract and view then next one or press key "2" to accept a contract.', 'center', scrn.ycenter+150, scrn.white);
        DrawFormattedText(instructions, 'If you have understood the instructions, press SPACE to continue', 'center', scrn.ycenter+200, scrn.white);
        
    end % end of phase statement 
    
    % copy the instructions window  and flip.
    Screen('CopyWindow',instructions,window,windrect, windrect);
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
    
    WaitSecs(1); % wait one sec before flipping to the block/trial/sequence information screen 
    
    %% ---------------------------------------
    % ADD THE TRIGGER INFORMATION (IF EEG = 1) 

    if EEG == 1

        % INIT COMMUNICATION WITH EXTERNAL DEVICES
        sp      = BioSemiSerialPort();
        set.sp  = sp;
        fprintf(' >>> OPENING USB TRIGGER LINK  <<<')

    end

    %% ---------------------------------------
    % START THE BLOCK & SEQUENCE/TRIAL LOOPS
    
    abort           = 0;                % when 1 subject can quit the experiment
    
    % UNPACK SETTINGS STRUCT
    ntrials         = set.totaltrials;  % total trials
    nb_blocks       = set.blocks;       % total blocks
    trialsPerBlock  = set.blocktrials;  % trials per block
    
    % INIT BLOCKS LOOP
    for iBlock = 1:nb_blocks
        
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

        if phase == 1 % run the phase 1 blocks 
            
            % UNPACK TRIALS STRUCT
            block_seq           = trials.sequence{iBlock};
            
            % display trial/sequence information window 
            Screen('OpenOffscreenWindow', window, windrect);
            Screen('TextSize', window, scrn.textsize);
            Screen('FillRect', window, scrn.grey ,windrect);
            DrawFormattedText(window, sprintf('Starting block %d',iBlock), 'center', scrn.ycenter-50, scrn.white);
            DrawFormattedText(window, 'Press SPACE to continue, or press ESC to quit', 'center', scrn.ycenter, scrn.white);
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
            if abort == 1 % exit 
                break;
            end
            
            set.iBlock      = iBlock;
            set.sequence    = block_seq;         % send the current sequence to the RUN function
            [set,logs]      = RunEconomic(set, scrn, logs);
            
        else % if phase == 2
            
            % UNPACK TRIALS STRUCT
            block_seq           = trials.sequence{iBlock};
            
            for trial = 1:trialsPerBlock
                
                set.iBlock = iBlock;
                set.thisTrial = trial;
                set.sequence    = block_seq{trial}; 
                
                % display trial/sequence information window 
                Screen('OpenOffscreenWindow', window, windrect);
                Screen('TextSize', window, scrn.textsize);
                Screen('FillRect', window, scrn.grey ,windrect);
                DrawFormattedText(window, sprintf('Starting sequence %d of block %d',trial, iBlock), 'center', scrn.ycenter-50, scrn.white);
                DrawFormattedText(window, 'Press SPACE to continue, or press ESC to quit', 'center', scrn.ycenter, scrn.white);
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
                if abort == 1 % exit 
                    break;
                end
                
                [set,logs]      = RunEconomic(set, scrn, logs); % run trials
                
                % UNPACK SET AND ADD THE TRIAL INFO TO THE "BLOCK"-LOG FILE 
                blocktrials(trial).session     = set.blocktrials.session;
                blocktrials(trial).block       = set.blocktrials.block;
                blocktrials(trial).trialnumber = set.blocktrials.trialnumber;
                blocktrials(trial).trialonset  = set.blocktrials.trialonset;
                blocktrials(trial).sequence    = set.blocktrials.sequence;
                blocktrials(trial).numsamples  = set.blocktrials.numsamples;
                blocktrials(trial).chosenitem  = set.blocktrials.chosenitem;

            end % End of trials loop
            
            % save trial info
            logs.blocktrials    = blocktrials;
            sub_log             = fullfile(logs.resultsfolder,sprintf(logs.blocktrialog,sub,taskName,iBlock,sess));
            save(sub_log,'logs');

        end % end of phase if statement
  
        % IF THIS IS THE LAST BLOCK BREAK FROM THE LOOP AND GO DIRECTLY TO
        % THE GOODBYE SCREEN
        if iBlock == nb_blocks
            break;
        end
        
        % if this is the end of the block 
        % allow subject to take a short break (if they want to)
        Screen('OpenOffscreenWindow', window, windrect);
        Screen('TextSize', window, scrn.textsize);
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
    end % end of blocks loop
    
    % THIS IS IT...
    % show thank you window
    Screen('OpenOffscreenWindow', window, windrect);
    Screen('TextSize', window, scrn.textsize);
    Screen('FillRect', window, scrn.grey ,windrect);
    DrawFormattedText(window, 'This is the end of the experiment. Thank you for your time', 'center', scrn.ycenter, scrn.white);
    Screen('Flip',window);
    WaitSecs(3);
    
    % clean up at the end of the experiment
    Screen('CloseAll');
    ShowCursor;
    Priority(0);
    fclose('all');
    
catch % ... catch last errors 
    Screen('CloseAll');
    ShowCursor;
    Priority(0);
    psychrethrow(psychlasterror);
end