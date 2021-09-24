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
taskNb          = 3; 

basedir         = pwd;

% add a new prompt to identify the phase of the task 
prompt          = {'Which phase is it?', 'Male or Female Stimuli?'};
dlgtitle        = 'Info window';
dims            = [1 30];
definput        = {'1', '1'}; % this is a default input (this should change)
taskphase       = inputdlg(prompt,dlgtitle,dims,definput);

phase           = str2num(taskphase{1});                                    % convert task phase to number
sess            = phase;
type            = str2num(taskphase{2});                                    % 1=female, 2=male, will be used to load the stimuli

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
    globalrect              = Screen('Rect', screenNumber);   
    
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
    scrn.globalrect         = globalrect;
    scrn.screenNumber       = screenNumber;
    
    %% ---------------------------------------
    % RUN A FEW IMPORTANT UTIL FUNCTIONS
    set                     = TaskSettings(taskNb, sess);                   % Define the first task-specific parameters
    
    set                     = DefineKeys(taskNb, set);                      % Define keys of the task
    
    set.type                = type;                                         % store type in settings struct to load the stimuli
    
    set                     = loaditems(set, wd);                           % read the excel file with the items (contracts)
   
    scrn.stimdeg            = set.stimsize_deg;
    scrn                    = screenSettings(scrn, taskNb);                 % Define screen setup
       
    [trials, set]           = CreateTrialList(set);                         % create trials, sequences, split in runs, etc..
    
    if phase == 2
        averaged                = GetAverage(taskName, wd, sub);            % get averaged ratimgs for each face
        set.items(:,2)          = averaged';                                % store the averaged ratings in the items list
    end
    %% ---------------------------------------
    % MAKE IMAGE TEXTURES FOR BOTH SIZES (NORMAL AND SMALL)
    % UNPACK STIMULI 
    data            = set.data;
    objects         = set.objects;
    
    % make a cell to store the image textures
    textures        = cell(1,objects);
    
    
    for i=1:objects
        
        % make textures
        textures{i} = Screen('MakeTexture', window, data(i).file); 
        
    end
    
    % UPDATE SETTIGS STRUCT
    set.textures = textures;
    
    % if phase is 2 make textures for the small images. The texture cell
    % and the smalltex cell will be desplayed together. When I try to
    % desplay all (main face and previous-small faces) using the texture
    % cell only I get an error. Apparently, this is how it works? 
    % the 
    if phase == 2
        
        % if this is the second phase, also include small textures
        smalltex        = cell(1,objects);       
        
        for j = 1:objects
            
            smalltex{j} = Screen('MakeTexture', window, data(j).file); 
        end
        
        set.smalltex = smalltex;
    end 
    
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
    
    if phase == 1 % show the instructions of phase 1 of the economic best choice task
        
        % display instructions for phase 1
        instructions = Screen('OpenOffscreenWindow', window, windrect);
        Screen('TextSize', instructions, scrn.textsize);
        Screen('FillRect', instructions, scrn.grey ,windrect);
        DrawFormattedText(instructions, 'This is the first phase of the experiment. You will be presented with images of faces', 'center', scrn.ycenter-200, scrn.white);
        DrawFormattedText(instructions, 'one-by-one at the centre of the screen. Your task is to carefully view each face and rate', 'center', scrn.ycenter-150, scrn.white);
        DrawFormattedText(instructions, '"how likely it would be for you to date that person in real life" on a scale of 1 to 9, where 1 means', 'center', scrn.ycenter-100, scrn.white);
        DrawFormattedText(instructions, '"I would never date that person" and 9 means "I would definitely date that person".','center', scrn.ycenter-50, scrn.white);
        DrawFormattedText(instructions, 'For you responses, press the keyboard keys 1 to 9 which correspond to your rating for each face.', 'center', scrn.ycenter, scrn.white);
        DrawFormattedText(instructions, 'Please rate the images of faces exactly as you would you in a real life scenario.', 'center', scrn.ycenter+50, scrn.white); 
        DrawFormattedText(instructions, 'If you have understood the instructions so far, press SPACE to continue', 'center', scrn.ycenter+100, scrn.white);
        
    else % if this is the second phase of the experiment 
        
        % display instructions for phase 2
        instructions = Screen('OpenOffscreenWindow', window, windrect);
        Screen('TextSize', instructions, scrn.textsize);
        Screen('FillRect', instructions, scrn.grey ,windrect);
        DrawFormattedText(instructions, 'This is the second phase of the experiment. On each trial/sequence of this phase,', 'center', scrn.ycenter-250, scrn.white);
        DrawFormattedText(instructions, 'you will be presented with up to 10 images of faces from the the previous phase, at the centre of the screen one-by-one.', 'center', scrn.ycenter-200, scrn.white);
        DrawFormattedText(instructions, 'Every time you are presented with a face, you may either "choose to date that person" or you may "reject it" and', 'center', scrn.ycenter-150, scrn.white);
        DrawFormattedText(instructions, 'view the next one. Please note that for each sequence you can choose only one face/date. If you reject','center', scrn.ycenter-100, scrn.white);
        DrawFormattedText(instructions, 'a face, you may not go back and select it. If, by the end of a sequence you have not chosen', 'center', scrn.ycenter-50, scrn.white);
        DrawFormattedText(instructions, 'a person to date, by default, the last face will be saved as your chosen date (for that sequence).', 'center', scrn.ycenter, scrn.white); 
        DrawFormattedText(instructions, 'Each face that you reject will be displayed at the bottom of the screen so that you have an idea of the', 'center', scrn.ycenter+50, scrn.white);
        DrawFormattedText(instructions, 'faces you rejected, and the number of faces left in the sequence. Press the keyboard key "1"', 'center', scrn.ycenter+100, scrn.white);
        DrawFormattedText(instructions, 'to accept a face/date or press the keyboard key "2" to reject the current face and view then next one.', 'center', scrn.ycenter+150, scrn.white);
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
        
        if phase == 1
            
            
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
            [set,logs]      = RunFace(set, scrn, logs);
            
        else % if this is phase 2
            
             % UNPACK TRIALS STRUCT
            block_seq           = trials.sequence{iBlock};
            
            for trial = 1:trialsPerBlock
                
                set.iBlock      = iBlock;
                set.thisTrial   = trial;
                set.sequence    = block_seq{trial};
                
                currentbalance  = set.balance;
                
                % display trial/sequence information window 
                Screen('OpenOffscreenWindow', window, windrect);
                Screen('TextSize', window, scrn.textsize);
                Screen('FillRect', window, scrn.grey ,windrect);
                DrawFormattedText(window, sprintf('Starting sequence %d of block %d',trial, iBlock), 'center', scrn.ycenter-50, scrn.white);
                DrawFormattedText(window, sprintf('Your current credit balance is %3.4f\n',currentbalance), 'center', scrn.ycenter, scrn.white);
                DrawFormattedText(window, 'Press SPACE to continue, or press ESC to quit', 'center', scrn.ycenter+50, scrn.white);
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
                
                [set,logs]      = RunFace(set, scrn, logs); % run trials
                
                % UNPACK SET AND ADD THE TRIAL INFO TO THE "BLOCK"-LOG FILE 
                blocktrials(trial).session      = set.blocktrials.session;
                blocktrials(trial).block        = set.blocktrials.block;
                blocktrials(trial).trialnumber  = set.blocktrials.trialnumber;
                blocktrials(trial).trialonset   = set.blocktrials.trialonset;
                blocktrials(trial).sequence     = set.blocktrials.sequence;
                blocktrials(trial).numsamples   = set.blocktrials.numsamples;
                blocktrials(trial).chosenitem   = set.blocktrials.chosenitem;
                blocktrials(trial).thisrate     = set.blocktrials.thisrate;
                blocktrials(trial).rank         = set.blocktrials.rank;
                blocktrials(trial).thisreward   = set.blocktrials.thisreward;
                blocktrials(trial).balance      = set.blocktrials.balance;
                
            end % End of trials loop
            
            % save trial info
            logs.blocktrials    = blocktrials;
            sub_log             = fullfile(logs.resultsfolder,sprintf(logs.blocktrialog,sub,taskName,iBlock,sess));
            save(sub_log,'logs');
            
        end % end of phase statement
        
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
        
    end % end of Blocks loop 
    
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
    
catch
    
    Screen('CloseAll');
    ShowCursor;
    Priority(0);
    psychrethrow(psychlasterror);
    
end