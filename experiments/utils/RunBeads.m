function [set, logs] = RunBeads(set, scrn, logs)

% THIS IS A SUBFUNCTION, PART OF THE "OPTIMAL STOPPING EXPERIMENTS". 

% it takes various stored information from other subfunctions to run
% one sequence

%% ---- Prepare the "global" information needed for all the tasks ---- %%

% UNPACK GLOBAL PARAMS FROM THE SETTINGS AND SCREEN STRUCTS
taskNb          = set.taskNb;       % number of task (needed to run the correct task)
thistrial       = set.thistrial;    % this is the current sequence/trial number 
thisblock       = set.thisblock;     % this is the current block number 

isi             = set.isi;          % interstimulus interval
jitter          = set.jitter; 
fixation        = set.fixation;     % draw fixation
EEG             = set.EEG;          % is it a EEG session? 

window          = scrn.window;       % main window
windrect        = scrn.windrect;
xcenter         = scrn.xcenter;
ycenter         = scrn.ycenter;
ifi             = scrn.ifi;          % frame duration
slack           = scrn.slack;        % slack is ifi/2 (important for timing)
white           = scrn.white;
black           = scrn.black;
grey            = scrn.grey;
globalrect      = scrn.globalrect;
fixsize         = scrn.fixationsize;
textfont        = scrn.textfont;
textsize        = scrn.textsize;
smalltext       = scrn.smalltext;

thisession      = logs.sess;
sub             = logs.sub;
taskname        = logs.task;
resfolder       = logs.resultsfolder;

% create fixation cross offscreen and paste later (faster)
fixationdisplay = Screen('OpenOffscreenWindow',window);
Screen('FillRect', fixationdisplay, grey);
Screen('TextFont',fixationdisplay, textfont);
Screen('TextSize',fixationdisplay, fixsize);
DrawFormattedText(fixationdisplay, fixation, 'center', ycenter, white);

%% ----- Run the beads task ------ %%

abort           = 0;
HideCursor;

% UNPACK BEADS RELATED PARAMETERS 
green           = scrn.green;
blue            = scrn.blue;
red             = scrn.red;

sequence        = set.sequence;     % this is the current sequence/trial 
thisurn         = set.urn;          % this is the current urn, contains the high prob colour beads 1=blue, 0=green
drawlen         = set.draws;        % length of each sequence (up to 10 draws)
penalty         = set.penalty;      % £0.25 loss for every draw
win             = set.win;          % £10 for every win!
balance         = set.balance;      % balance starts as zero and updates from correct/incorrect (loss) responses & draw choices

bead_dur        = set.bead_dur;     % duration of each bead in sec
response        = set.response;     % response duration in sec
fix_dur         = set.fix_dur;      % fixation duration
dfeedback       = set.feed_dur;     % duration of feedback window in sec

bluekey         = set.code1;       % keycode for blue option
greenkey        = set.code2;       % keycode for green option
drawkey         = set.code3;       % keycode for draw-again option
esckey          = set.code21;      % keycode for aborting the experiment 

anchors         = set.anchors;
scalelength     = set.scalelength;
maxtime         = set.maxtime;
mousebutton     = set.mousebutton;
leftclick       = set.leftclick;
rightclick      = set.rightclick;
midclick        = set.midclick;
horzline        = set.horzline;
width           = set.width;
scalepos        = set.scalepos;
startpos        = set.startpos;
line            = set.line;
scalerange      = set.scalerange;
scalecolour     = white;
slidercolour    = black;

% calculate coordinates of scale line and text bounds
x               = globalrect(3)/2;

% ADD A FEW MORE VARIABLES
blocktrials     = []; % store trial info
draws           = []; % store info specifically about the sequence/draws [1-10]
previous        = []; % store the previous draws here to show at the bottom of the screen
colours         = []; % store the previous draw-colours here to show at the bottom of the screen
draw_count      = 0;  % drawing counter
accuracy        = nan;

% UNPACK EEG TRIGGERS
if EEG == 1 
    
    sp          = set.sp;
    trigger1    = set.trigger1;
    trigger2    = set.trigger2;
    trigger3    = set.trigger3;
    trigger4    = set.trigger4;
    trigger5    = set.trigger5;
    
    trigger6    = set.trigger6;
    trigger7    = set.trigger7;
    trigger8    = set.trigger8;
    trigger9    = set.trigger9;
    trigger10   = set.trigger10;
    trigger11   = set.trigger11;
    trigger12   = set.trigger12;
    trigger13   = set.trigger13;
    trigger14   = set.trigger14;
    trigger15   = set.trigger15;
    trigger16   = set.trigger16;
    trigger17   = set.trigger17;
    
    trigger102  = set.trigger102; % start of sequence
    trigger103  = set.trigger103; % end of sequence
    
end
    
% SET ALL RELEVANT WINDOWS (BEAD WINDOWS, RESPONSE WINDOW, CONFIDENCE RATING & FEEDBACK WINDOWS)

% BLUE WINDOW
blue_window = Screen('OpenOffscreenWindow',window);
Screen('TextSize', blue_window, textsize);
Screen('FillRect', blue_window, grey ,windrect);
DrawFormattedText(blue_window, 'BLUE', 'center', ycenter, blue);

% GREEN WINDOW
green_window = Screen('OpenOffscreenWindow',window);
Screen('TextSize', green_window, textsize);
Screen('FillRect', green_window, grey ,windrect);
DrawFormattedText(green_window, 'GREEN', 'center', ycenter, green);

% RESPONSE WINDOW (with options blue, green and draw-again)
response_window1 = Screen('OpenOffscreenWindow',window);
Screen('TextSize', response_window1, textsize);
Screen('FillRect', response_window1, grey ,windrect);
DrawFormattedText(response_window1, 'BLUE?', 'center', ycenter-50, blue);
DrawFormattedText(response_window1, 'GREEN?', 'center', ycenter, green);
DrawFormattedText(response_window1, 'DRAW?', 'center', ycenter+50, white);

% RESPONSE WINDOW (with options blue, green)
response_window2 = Screen('OpenOffscreenWindow',window);
Screen('TextSize', response_window2, textsize);
Screen('FillRect', response_window2, grey ,windrect);
DrawFormattedText(response_window2, 'That was your last draw:', 'center', ycenter-50, white);
DrawFormattedText(response_window2, 'BLUE?', 'center', ycenter, blue);
DrawFormattedText(response_window2, 'GREEN?', 'center', ycenter+50, green);

% FEEDBACK WINDOW (you win!)
feedback_window1 = Screen('OpenOffscreenWindow',window);
Screen('TextSize', feedback_window1, textsize);
Screen('FillRect', feedback_window1, grey ,windrect);
DrawFormattedText(feedback_window1, 'Well Done! :)', 'center', ycenter, green);

% FEEDBACK WINDOW (you lose!)
feedback_window2 = Screen('OpenOffscreenWindow',window);
Screen('TextSize', feedback_window2, textsize);
Screen('FillRect', feedback_window2, grey ,windrect);
DrawFormattedText(feedback_window2, 'Oh No! :(', 'center', ycenter, red);

% FEEDBACK WINDOW (wrong answer!)
feedback_window3 = Screen('OpenOffscreenWindow',window);
Screen('TextSize', feedback_window3, textsize);
Screen('FillRect', feedback_window3, grey ,windrect);
DrawFormattedText(feedback_window3, 'Oh No! :(', 'center', ycenter-50, red);
DrawFormattedText(feedback_window3, 'You you are not allowed to draw again', 'center', ycenter, red);

% FEEDBACK WINDOW (didn't respond!)
feedback_window4 = Screen('OpenOffscreenWindow',window);
Screen('TextSize', feedback_window4, textsize);
Screen('FillRect', feedback_window4, grey ,windrect);
DrawFormattedText(feedback_window4, 'Oh No! :(', 'center', ycenter-50, red);
DrawFormattedText(feedback_window4, 'You did not respond', 'center', ycenter, red);

% ADD TEXTBOUNDS -- WILL BE USED TO CREATE THE SLIDER 
textbounds = [Screen('TextBounds', window, sprintf(anchors{1})); Screen('TextBounds', window, sprintf(anchors{3}))];
         
% CONFIDENCE RATING WINDOW
confidence_window = Screen('OpenOffscreenWindow',window);
Screen('TextSize', confidence_window, textsize);
Screen('FillRect', confidence_window, grey ,windrect);
DrawFormattedText(confidence_window, 'On a scale of 0 to 100, how confident are you for your choice?', 'center', ycenter-150, white);
DrawFormattedText(confidence_window, '0 = Not Confiddent', 'center', ycenter-50, white);
DrawFormattedText(confidence_window, '100 = Very Confident', 'center', ycenter, white);

% Left, middle and right anchors
DrawFormattedText(confidence_window, anchors{1}, leftclick(1, 1) - textbounds(1, 3)/2,  windrect(4)*scalepos+40, [],[],[],[],[],[],[]); % Left point
DrawFormattedText(confidence_window, anchors{2}, 'center',  windrect(4)*scalepos+40, [],[],[],[],[],[],[]); % Middle point
DrawFormattedText(confidence_window, anchors{3}, rightclick(1, 1) - textbounds(2, 3)/2,  windrect(4)*scalepos+40, [],[],[],[],[],[],[]); % Right point

% Drawing the scale
Screen('DrawLine', confidence_window, scalecolour, midclick(1), midclick(2), midclick(3), midclick(4), width);         % Mid tick
Screen('DrawLine', confidence_window, scalecolour, leftclick(1), leftclick(2), leftclick(3), leftclick(4), width);     % Left tick
Screen('DrawLine', confidence_window, scalecolour, rightclick(1), rightclick(2), rightclick(3), rightclick(4), width); % Right tick
Screen('DrawLine', confidence_window, scalecolour, horzline(1), horzline(2), horzline(3), horzline(4), width);     % Horizontal line

%%
% START THE TRIAL WITH A FIXATION CROSS
Screen('CopyWindow', fixationdisplay,window, windrect, windrect)
fliptime    = Screen('Flip', window); % flip fixation window
trialstart  = fliptime;

% send sequence start trigger
if EEG == 1 
    sp.sendTrigger(trigger102)
end

% object offset
object_offset   = trialstart + fix_dur + isi + randperm(jitter*1000,1)/1000 - ifi;

% 1. BEGIN DRAWING 
for thisdraw = 1:drawlen
    
    % allow subject to abort the experiment 
    [keyisdown,~,keycode] = KbCheck;
    if keyisdown && keycode(esckey)
        abort = 1;
        break;
    end
    
    xcenters        = xcenter - 250; % start adding the darws prices at xcenter - 600
    xs              = xcenters;
    
    % STATEMENT 1. Is it a high or a low probability draw?
    if sequence(thisdraw) == 1  % high prob draw
        
        % STATEMENT 2. Is it a blue or a green bead?
        if thisurn == 1         % blue urn and blue bead
            % show blue window
            Screen('CopyWindow', blue_window, window, windrect, windrect)
            if EEG == 1
                stimtrigger         = trigger1; % assign trigger 1 as stimulus triger
            end
            previous_colour     = blue;
            previous_bead       = 'blue';
            
        else % green urn and green beed
            % show green window
            Screen('CopyWindow', green_window, window, windrect, windrect)
            if EEG == 1
                stimtrigger         = trigger3;
            end
            previous_colour     = green;
            previous_bead       = 'green';
         
        end
            
    elseif sequence(thisdraw) == 2 % low prob draw
        
        if thisurn == 1 % blue urn and green bead
            % show green window
            Screen('CopyWindow', green_window, window, windrect, windrect)
            if EEG == 1
                stimtrigger         = trigger2;
            end
            previous_colour     = green;
            previous_bead       = 'green';
            
        else % green urn and blue bead
            % show blue window
            Screen('CopyWindow', blue_window, window, windrect, windrect)
            if EEG == 1
                stimtrigger         = trigger4;
            end
            previous_colour     = blue;
            previous_bead       = 'blue';
            
        end 
    end % end of statment 1 if 
    
    % if this is not the first draw, show previous draw at the bottom of
    % the screen
    if thisdraw > 1
        
        previous_len        = length(previous); % how many previous prices?
        
        for i = 1:previous_len
            
            % show the previous prices at the bottom of the screen
            Screen('TextSize', window, smalltext);
            DrawFormattedText(window, 'Previous Draws', xcenter - 250, ycenter+250, white); 
            DrawFormattedText(window, previous{i}, xcenters, ycenter+350, colours{i}); 

            % update xcenters, so that previous prices are not
            % displayed on top of each other
            xcenters = xcenters + 80;
            
        end
    end
    
    bead_onset     = Screen('Flip', window, object_offset - slack); 
    bead_offset    = bead_onset + bead_dur - ifi;                          % bead on for 0.5 ms
    
    % send sequence start trigger
    if EEG == 1
        sp.sendTrigger(stimtrigger) % blue urn -- high prob blue bead trigger
    end
    
    % 2. BRING FIXATION BACK ON
    Screen('CopyWindow', fixationdisplay,window, windrect, windrect)
  
    % if this is not the first draw, show previous draw at the bottom of
    % the screen (during fixation also
    if thisdraw > 1
        
        previous_len        = length(previous); % how many previous prices?
        
        for i = 1:previous_len
            
            % show the previous prices at the bottom of the screen
            Screen('TextSize', window, smalltext);
            DrawFormattedText(window, 'Previous Draws', xcenter - 250, ycenter+250, white); 
            DrawFormattedText(window, previous{i}, xs, ycenter+350, colours{i}); 

            % update xcenters, so that previous prices are not
            % displayed on top of each other
            xs = xs + 80;
            
        end
    end
    
    fixation_onset      = Screen('Flip', window, bead_offset - slack);        % fixation on, prepare for response prompt  
    object_offset       = fixation_onset + fix_dur - ifi;     % add jitter here?
    
    % 3. SHOW RESPONSE PROMPT 
    if thisdraw < drawlen % show response window 1
        Screen('CopyWindow', response_window1, window, windrect, windrect)
        
    else % show response window 2
        Screen('CopyWindow', response_window2, window, windrect, windrect)
      
    end % end of resp prompt if statement
    
    prompt_onset    = Screen('Flip', window, object_offset - slack); 
    
    % send response prompt trigger
    if EEG == 1 
        sp.sendTrigger(trigger5)
    end
    
    % 3. INIT RESPONSE: initiate on each draw
    rt                  = NaN;
    answer              = NaN;
    resp_input          = 0;
    responseTrigNotSent = 1;
    
    while resp_input == 0 && (GetSecs - prompt_onset) < response 
        
        [keyisdown, secs, keycode] = KbCheck;
        pressedKeys = find(keycode);
        
        % check the response key
        if isempty(pressedKeys) % if something weird happens or subject doesn't press any of the correct keys
            resp_input  = 0; 
            rt          = nan;
            answer      = nan;
            respmade    = GetSecs;
            
        elseif ~isempty(pressedKeys) % if subject pressed a valid key
            
            if keycode(1,bluekey) % if subject chose the blue urn 
                resp_input  = bluekey;
                rt          = secs - prompt_onset;
                answer      = 1; % blue urn 
                respmade    = secs;
                
                % send blue choice trigger
                if EEG == 1 && responseTrigNotSent==1
                    sp.sendTrigger(trigger6);
                    responseTrigNotSent=0;
                end
                
            elseif keycode(1,greenkey) % if subject chose the green urn 
                resp_input  = greenkey;
                rt          = secs - prompt_onset;
                answer      = 2; % green urn
                respmade    = secs;
                
                % send green choice trigger
                if EEG == 1 && responseTrigNotSent==1
                    sp.sendTrigger(trigger7);
                    responseTrigNotSent=0;
                end
                
            elseif keycode(1,drawkey) % if subject chose to draw again
                resp_input  = drawkey;
                rt          = secs - prompt_onset;
                answer      = 3; % draw-again 
                respmade    = secs;
                draw_count  = draw_count + 1;
                
                % send draw again trigger
                if EEG == 1 && responseTrigNotSent==1
                    sp.sendTrigger(trigger8);
                    responseTrigNotSent=0;
                end
            end %  
        end % if responded statement
    end % end of response while loop 
    
    prompt_offset   = respmade + isi - ifi;                                     % response prompt self paced or on for 2500 ms
    
   % 4. BRING BACK FIXATION. NEXT STEP:
   % A) if answer = 3, draw again,
   % B) if answer = 1 or answer = 2, ask subject to rate the confidence of their choice,
   % check if response was correct or incorrect, give appropreate feedback, break 
   % C) if draw_count == 10 and subject chose to draw again collect it as
   % wrong answer 
   
   Screen('CopyWindow', fixationdisplay,window, windrect, windrect)
   fixation_onset   = Screen('Flip', window, prompt_offset - slack);                % fixation on, prepare for next draw or for feedback      
  
   object_offset    = fixation_onset + isi + randperm(jitter*1000,1)/1000 - ifi;     % add jitter here?
   
   if answer ~= 3
       
       if any(isnan(answer)) % if subject doesn't respond, just move to the next draw
           
           % show feedback screen 4 and move to next trial
           Screen('CopyWindow', feedback_window4, window, windrect, windrect)
           feedback_onset = Screen('Flip', window, object_offset - slack);    % rating window is on
           
           if EEG == 1
               sp.sendTrigger(trigger17);
           end
           
           % feedback is on
           object_offset = feedback_onset + dfeedback - ifi;
           
            % BRING FIXATION BACK ON
           Screen('CopyWindow', fixationdisplay,window, windrect, windrect)
           object_onset     = Screen('Flip', window, object_offset - slack); 

           object_offset    = object_onset + isi + randperm(jitter*1000,1)/1000 - ifi;% add jitter here?
           
       else
           
           % THIS IS THE PART WHERE WE DRAW THE SLIDER AND SHOW THE
           % CONFIDENCE SCREEN
           
           % First define the starting point of the slider
           % calculate coordinates of scale line and text bounds
           if strcmp(startpos, 'left')
               x = globalrect(3)*(1-scalelength);
           elseif strcmp(startpos, 'center')
               x = globalrect(3)/2;
           elseif strcmp(startpos, 'right')
               x = globalrect(3)*scalelength;
           end
           
           % NOW DARW THE SCALE WITHOUT THE SLIDER AND ALLOW THE SUBJECT TO
           % CLICK ONCE TO INIT THE SLIDER
           
           % initialise the mouse
           SetMouse(round(x), round(windrect(4)*scalepos), window, 1)
           
           t0                         = GetSecs;
           initresp                   = 0;
           responseTrigNotSent        = 1;
           
           while initresp == 0 && (GetSecs - fixation_onset) < maxtime
               
               [x,~,buttons,~,~,~] = GetMouse(window, 1);
               
               % Stop at upper and lower bound
               if x > windrect(3)*scalelength
                   x = windrect(3)*scalelength;
               elseif x < windrect(3)*(1-scalelength)
                   x = windrect(3)*(1-scalelength);
               end
               
               % draw the question and slider
               Screen('CopyWindow', confidence_window, window, windrect, windrect)
               rating_onset = Screen('Flip', window, object_offset - slack);    % rating window is on
               
                % wait for second response 
                secs = GetSecs;
                if buttons(mousebutton) == 1
                    initresp = 1;
                end

                % Abort if answer takes too long
                if secs - t0 > maxtime 
                    break
                end
           end % end of first mouse while loop
           
           KbReleaseWait; %Keyboard
           WaitSecs(0.1) % % delay to prevent fast mouse clicks mix 
           
           % NOW that the saubject has clicked the mouse once or enough
           % time has passed, aske them to give their rate by dragging the
           % mouse 
           t1               = GetSecs;
           secondresp       = 0;
           
           while secondresp == 0 && (GetSecs - rating_onset) < maxtime
               
               [x,~,buttons,~,~,~] = GetMouse(window, 1);
               
               % Stop at upper and lower bound
               if x > windrect(3)*scalelength
                   x = windrect(3)*scalelength;
               elseif x < windrect(3)*(1-scalelength)
                   x = windrect(3)*(1-scalelength);
               end
               
               % draw the question and slider
               Screen('TextSize', window, textsize);
               Screen('FillRect', window, grey ,windrect);
               DrawFormattedText(window, 'On a scale of 0 to 100, how confident are you for your choice?', 'center', ycenter-150, white);
               DrawFormattedText(window, '0 = Not Confiddent', 'center', ycenter-50, white);
               DrawFormattedText(window, '100 = Very Confident', 'center', ycenter, white);

               % Left, middle and right anchors
               DrawFormattedText(window, anchors{1}, leftclick(1, 1) - textbounds(1, 3)/2,  windrect(4)*scalepos+40, [],[],[],[],[],[],[]); % Left point
               DrawFormattedText(window, anchors{2}, 'center',  windrect(4)*scalepos+40, [],[],[],[],[],[],[]); % Middle point
               DrawFormattedText(window, anchors{3}, rightclick(1, 1) - textbounds(2, 3)/2,  windrect(4)*scalepos+40, [],[],[],[],[],[],[]); % Right point

                % Drawing the scale
               Screen('DrawLine', window, scalecolour, midclick(1), midclick(2), midclick(3), midclick(4), width);         % Mid tick
               Screen('DrawLine', window, scalecolour, leftclick(1), leftclick(2), leftclick(3), leftclick(4), width);     % Left tick
               Screen('DrawLine', window, scalecolour, rightclick(1), rightclick(2), rightclick(3), rightclick(4), width); % Right tick
               Screen('DrawLine', window, scalecolour, horzline(1), horzline(2), horzline(3), horzline(4), width);     % Horizontal line

               % The slider
               Screen('DrawLine', window, slidercolour, x, windrect(4)*scalepos - line, x, windrect(4)*scalepos  + line, width);

               position = round((x)-min(scalerange));                       % Calculates the deviation from 0. 
               position = (position/(max(scalerange)-min(scalerange)))*100; % Converts the value to percentage

               DrawFormattedText(window, num2str(round(position)), 'center', windrect(4)*(scalepos - 0.05), white); 
               rating_onset = Screen('Flip', window);    % rating window is on
               
               % send confidence screen trigger
               if EEG == 1 
                   sp.sendTrigger(trigger9)
               end
               
               secs = GetSecs;
               if buttons(mousebutton) == 1
                   secondresp = 1;
               end
               rate_rt = (secs - t1);

               % sent rate (response) trigger
               if EEG == 1 && responseTrigNotSent==1
                   
                   if round(position) <= 25 
                       sp.sendTrigger(trigger10);
                       
                   elseif round(position) > 25 && round(position) <= 50 
                       sp.sendTrigger(trigger11);
                       
                   elseif round(position) > 50 && round(position) <= 75
                       sp.sendTrigger(trigger12);
                       
                   else 
                       sp.sendTrigger(trigger13);
                   end 
                   responseTrigNotSent=0;
               end
               
           end % end of second mouse while loop
           
           object_offset   = secs + isi - ifi;
          
           KbReleaseWait; %Keyboard
           
           % ADD A FIXATION BEFORE FEEDBACK
           Screen('CopyWindow', fixationdisplay, window, windrect, windrect)
           object_onset     = Screen('Flip', window, object_offset - slack); 

           object_offset    = object_onset + fix_dur - ifi;
           

           if answer == 1 % if subject chose blue urn
               if thisurn == 1                                                      % if thisurn is in fact a blue urn

                   accuracy         = 1;                                            % update the accuracy var

                   Screen('CopyWindow', feedback_window1,window, windrect, windrect)% show "you win" feedback window
                   feedback_onset = Screen('Flip', window, object_offset - slack);  % feedback window on 

                   % send you win trigger
                    if EEG == 1
                        sp.sendTrigger(trigger14);
                    end

               else % if "thisurn" is a green urn 
                   accuracy         = 0;

                   Screen('CopyWindow', feedback_window2,window, windrect, windrect) % show "you lose" feedback window
                   feedback_onset = Screen('Flip', window, object_offset - slack);   % feedback window on 

                   % send you lose trigger
                    if EEG == 1
                        sp.sendTrigger(trigger15);
                    end
               end
               object_offset = feedback_onset + dfeedback - ifi;
               
               % show fixation between ending the trial/sequence
               Screen('CopyWindow', fixationdisplay,window, windrect, windrect)
               object_onset     = Screen('Flip', window, object_offset - slack); 

               object_offset    = object_onset +  isi + randperm(jitter*1000,1)/1000 - ifi; % add jitter here?
               break;

           elseif answer == 2 % if subject chose green urn 

               if thisurn == 1 % if thisurn is a blue urn

                   accuracy         = 0;

                   Screen('CopyWindow', feedback_window2,window, windrect, windrect) % show "you lose" feedback window
                   feedback_onset = Screen('Flip', window, object_offset - slack);   % feedback window on 
                   % send you lose trigger
                    if EEG == 1
                        sp.sendTrigger(trigger15);
                    end

               else % if "thisurn" is a green urn 

                   Screen('CopyWindow', feedback_window1,window, windrect, windrect) % show "you win" feedback window
                   accuracy = 1;                                                     % update the accuracy var
                   feedback_onset = Screen('Flip', window, object_offset - slack);   % feedback window on 
                   % send you lose trigger
                   if EEG == 1
                       sp.sendTrigger(trigger14);
                   end
               end
               object_offset = feedback_onset + dfeedback - ifi;
               
               % show fixation between ending the trial/sequence
               Screen('CopyWindow', fixationdisplay,window, windrect, windrect)
               object_onset     = Screen('Flip', window, object_offset - slack); 

               object_offset    = object_onset +  isi + randperm(jitter*1000,1)/1000 - ifi; % add jitter here?
               break

           end
       end
       
   elseif answer == 3
       
       if draw_count == drawlen % if this is the last (10th) draw ans subject requested another draw
           
           accuracy     = 0;
           % DISPLAY THE 3RD FEEDBACK SCREEN
           Screen('CopyWindow', feedback_window3,window, windrect, windrect)
           feedback_onset   = Screen('Flip', window, object_offset - slack); % feedback window on 
           % send you lose (out of draws) trigger
           if EEG == 1
               sp.sendTrigger(trigger16);
           end
           
           % update answer, given that participant pressed 3 (for draw 10), answer should be 0 (this is an incorrect trial)
           answer           = 0;
           
           object_offset    = feedback_onset + dfeedback - ifi;
           
           % BRING FIXATION BACK ON BEFORE ENDING THE TRIAL
           Screen('CopyWindow', fixationdisplay,window, windrect, windrect)
           object_onset     = Screen('Flip', window, object_offset - slack); 
           
           object_offset    = object_onset +  isi + randperm(jitter*1000,1)/1000 - ifi; % add jitter here?
           break
           
       else % if this is not the last draw and subject wants to draw again
           
           Screen('CopyWindow', fixationdisplay,window, windrect, windrect)
           object_onset     = Screen('Flip', window, object_offset - slack); 
           
           object_offset    = object_onset +  isi + randperm(jitter*1000,1)/1000 - ifi; % add jitter here?
           
       end
   end
   
   % store the current price to show at the bottom of the screen
   % (during the next samples)
   previous{thisdraw}          = previous_bead;
   colours{thisdraw}           = previous_colour;
           
   % add sequence/draws related info here 
   draws(thisdraw).session      = thisession;
   draws(thisdraw).block        = thisblock;
   draws(thisdraw).trialnumber  = thistrial;
   draws(thisdraw).trialonset   = trialstart;
   draws(thisdraw).thisdraw     = thisdraw;
   draws(thisdraw).rt           = rt;
   
   if abort; fclose('all');break; end 

end % end of draw for loop

% Do we send a trigger at this point? something like:
if EEG == 1 
    sp.sendTrigger(trigger103)

end

% UPDATE BALANCE 
% 1. First subtract the penalty (0.25p for every draw)
balance = balance - (penalty * draw_count);

% 2. Now add winnings or subtract losses
if accuracy == 1
    balance = balance + win;
  
elseif accuracy == 0 
    balance = balance - win;
    
end

% update the current balance
set.balance         = balance;

% update the position/rate for this trial (if the subject didn't give a
% rate then this should be nan
if answer == 1 || answer == 2
    thisrate        = position;
    
else 
    position        = nan;
    thisrate        = position;
    rate_rt         = nan;
end
   
% add trial related info here 
blocktrials.session      = thisession;
blocktrials.block        = thisblock;
blocktrials.trialnumber  = thistrial;
blocktrials.trialonset   = trialstart;
blocktrials.urntype      = thisurn;
blocktrials.sequence     = sequence(1,:);
blocktrials.draws        = draw_count;
blocktrials.response     = answer;
blocktrials.accuracy     = accuracy;
blocktrials.balance      = balance;
blocktrials.thisrate     = thisrate;
blocktrials.ratert       = rate_rt;

WaitSecs(2); %  1= force wait for actual pulse; 0=return this many ms after pulse

% store draws and trials info in the logs mat file 
logs.trialstart         = trialstart;
logs.draws              = draws;
set.blocktrials         = blocktrials;

sub_drawslog             = fullfile(resfolder,sprintf(logs.drawslog,sub,taskname,thisblock,thistrial,thisession));
save(sub_drawslog,'logs');


% end of run/sequence 


end