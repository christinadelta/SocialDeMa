function [set, logs] = RunTrials(set, scrn, logs, keys)

% THIS IS A SUBFUNCTION, PART OF THE "OPTIMAL STOPPING EXPERIMENTS". 

% it takes various stored information from other subfunctions to run
% one sequence

%% ---- Prepare the "global" information needed for all the tasks ---- %%

% UNPACK GLOBAL PARAMS FROM THE SETTINGS AND SCREEN STRUCTS
taskNb          = set.taskNb;       % number of task (needed to run the correct task)
total_trials    = set.trials;       % total number of trials

isi             = set.isi;          % interstimulus interval
jitter          = set.jitter; 
fixation        = set.fixation;     % draw fixation

window          = scrn.window;       % main window
windrect        = scrn.windrect;
xcenter         = scrn.xcenter;
ycenter         = scrn.ycenter;
ifi             = scrn.ifi;          % frame duration
slack           = scrn.slack;        % slack is ifi/2 (important for timing)
white           = scrn.white;
grey            = scrn.grey;
black           = scrn.black;
fixsize         = scrn.fixationsize;
textfont        = scrn.textfont;
textsize        = scrn.textsize;


% keys            = DefineKeys(taskNb); % run the keys function


% create fixation cross offscreen and paste later (faster)
fixationdisplay = Screen('OpenOffscreenWindow',window);
Screen('FillRect', fixationdisplay, grey);
Screen('TextFont',fixationdisplay, textfont);
Screen('TextSize',fixationdisplay, fixsize);
DrawFormattedText(fixationdisplay, fixation, xcenter, ycenter, white);

%% ----- Run the correct task ------ %%

% UNPACK BEADS RELATED PARAMETERS 
green       = scrn.green;
blue        = scrn.blue;
red         = scrn.red;

sequence    = set.sequence; % this is the current sequence/trial 
thistrial   = set.thistrial; % this is the current sequence/trial number 
thisurn     = set.urn; % this is the current urn, contains the high prob colour beads 1=blue, 0=green
drawlen     = set.draws; % length of each sequence (up to 10 draws)

bead_dur    = set.bead_dur; % duration of each bead in sec
response    = set.response; % response duration in sec
dfeedback   = set.feed_dur; % duration of feedback window in sec

bluekey     = keys.code1; % keycode for blue option
greenkey    = keys.code2; % keycode for green option
drawkey     = keys.code3; % keycode for draw-again option

% SET ALL RELEVANT WINDOWS (BEAD WINDOWS, RESPONSE WINDOWS)

% BLUE WINDOW
blue_window = Screen('OpenOffscreenWindow',window);
Screen('TextSize', blue_window, textsize);
Screen('FillRect', blue_window, grey ,windrect);
DrawFormattedText(blue_window, 'BLUE', 'center', 'center', blue);

% GREEN WINDOW
green_window = Screen('OpenOffscreenWindow',window);
Screen('TextSize', green_window, textsize);
Screen('FillRect', green_window, grey ,windrect);
DrawFormattedText(green_window, 'GREEN', 'center', 'center', green);

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
DrawFormattedText(response_window2, 'That was your last draw:', 'center', ycenter-50, blue);
DrawFormattedText(response_window2, 'BLUE?', 'center', ycenter, blue);
DrawFormattedText(response_window2, 'GREEN?', 'center', ycenter+50, green);

% FEEDBACK WINDOW (you win!)
feedback_window1 = Screen('OpenOffscreenWindow',window);
Screen('TextSize', feedback_window1, textsize);
Screen('FillRect', feedback_window1, grey ,windrect);
DrawFormattedText(feedback_window1, 'Well Done! :)', 'center', ycenter-50, green);
DrawFormattedText(feedback_window1, 'You win £10', 'center', ycenter, green);

% FEEDBACK WINDOW (you lose!)
feedback_window2 = Screen('OpenOffscreenWindow',window);
Screen('TextSize', feedback_window2, textsize);
Screen('FillRect', feedback_window2, grey ,windrect);
DrawFormattedText(feedback_window2, 'Oh No! :(', 'center', ycenter-50, red);
DrawFormattedText(feedback_window2, 'You lose £10', 'center', ycenter, red);

% FEEDBACK WINDOW (wrong answer!)
feedback_window3 = Screen('OpenOffscreenWindow',window);
Screen('TextSize', feedback_window3, textsize);
Screen('FillRect', feedback_window3, grey ,windrect);
DrawFormattedText(feedback_window3, 'Wrong! :(', 'center', ycenter-50, red);
DrawFormattedText(feedback_window3, 'You you are not allowd to draw again', 'center', ycenter, red);

trials      = [];   % store trial info
draw_count  = 0;    % drawing counter

% START THE TRIAL WITH A FIXATION CROSS
Screen('CopyWindow', fixationdisplay,window, windrect, windrect)
fliptime    = Screen('Flip', window); % flip fixation window
trialstart  = fliptime;

% object offset
objectoff   = trialstart + isi + randperm(jitter*1000,1)/1000 - ifi;

% 1. BEGIN DRAWING 
for thisdraw = 1:drawlen
    
    % STATEMENT 1. Is it a high or a low probability draw?
    if sequence(thisdraw) == 1 % high prob draw
        
        % STATEMENT 2. Is it a blue or a green bead?
        if thisurn == 1 % blue urn and blue bead
            % show blue window
            Screen('CopyWindow', blue_window, window, windrect, windrect)
            
        else % green urn and green beed
            % show green window
            Screen('CopyWindow', green_window, window, windrect, windrect)
        end
            
    elseif sequence(thisdraw) == 2 % low prob draw
        
        if thisurn == 1 % blue urn and green bead
            % show green window
            Screen('CopyWindow', green_window, window, windrect, windrect)
            
        else % green urn and blue bead
            % show blue window
            Screen('CopyWindow', blue_window, window, windrect, windrect)
            
        end 
    end % end of statment 1 if 
    
    beadon      = Screen('Flip', window, objectoff - slack);            % here the current bead is fliped
    objectstart = beadon;
    
    beadstart   = objectstart - trialstart;                             % timestamp -- start of the first bead
    beadoff     = beadon + bead_dur - ifi;                            % bead on for 1000 ms
    
    % 2. SHOW RESPONSE PROMPT 
    if thisdraw < drawlen % show response window 1
        Screen('CopyWindow', response_window1, window, windrect, windrect)
        
    else % show response window 2
        Screen('CopyWindow', response_window2, window, windrect, windrect)
      
    end % end of resp prompt if statement
    
    prompton    = Screen('Flip', window, beadoff - slack); 
    
    fprintf('bead was on for %3.4f\n', prompton - objectstart);        % time interval from the first flip until the next
    
    
    % 3. INIT RESPONSE: initiate on each draw
    rt          = NaN;
    answer      = NaN;
    resp_input  = 0;
    
    while resp_input == 0 && (GetSecs - prompton) < response 
        
        [keyisdown, secs, keycode] = KbCheck;
        pressedKeys = find(keycode);
        
        % check the response key
        if isempty(pressedKeys) % if something weird happens or subject doesn't press any of the correct keys
            resp_input  = 0; 
            rt          = nan;
            answer      = nan;
            respmade    = GetSecs;
            
        elseif ~isempty(pressedKeys) % if subject pressed a valid key
            
            if keycode(1,bluekey) % subject if subject chose the green urn 
                resp_input  = bluekey;
                rt          = secs - prompton;
                answer      = 1; % green urn 
                respmade    = secs;
                
            elseif keycode(1,greenkey) % subject pressed the inanimate key
                resp_input  = greenkey;
                rt          = secs - prompton;
                answer      = 2; % blue
                respmade    = secs;
                
            elseif keycode(1,drawkey) % subject pressed the inanimate key
                resp_input  = drawkey;
                rt          = secs - prompton;
                answer      = 3; % draw-again 
                respmade    = secs;
                draw_count  = draw_count + 1; 
            end %  % need to add escape key here 
        end % if responded statement
    end % end of response while loop 
    
    promptoff   = prompton + response - ifi;                            % response prompt on for 2500 ms
    set.answer  = answer;
    
   % 3. BRING BACK FIXATION. NEXT STEP:
   % A) if answer = 3, draw again,  
   % B) if answer = 1 or answer = 2, check if response was correct or
   % incorrect, give appropreate feedback, break 
   % C) if draw_count == 10 and subject chose to draw again c0llect it as
   % wrong answer and breat 
   
   Screen('CopyWindow', fixationdisplay,window, windrect, windrect)
   fixon        = Screen('Flip', window, promptoff - slack);                % fixation on, prepare for next draw or for feedback      
   
   fprintf('prompt was on for %3.4f\n', fixon - prompton);                  % time interval from the the flip of the prompt until the next
   
   objectoff    = fixon + isi + randperm(jitter*1000,1)/1000 - ifi;         % add jitter here?
   
   % fprintf('fixation was on for %3.4f\n', objectoff); 
   
   if answer ~= 3
       
       if answer == 1 % if subject chose blue urn
           if thisurn == 1 % if thisurn is in fact a blue urn
           
               Screen('CopyWindow', feedback_window1,window, windrect, windrect) % show "you win" feedback window
               

           else % if "thisurn" is a green urn 

               Screen('CopyWindow', feedback_window2,window, windrect, windrect) % show "you lose" feedback window

           end
           
           feedbackon = Screen('Flip', window, objectoff - slack);      % feedback window on 
           
           fprintf('fixation was on for %3.4f\n', feedbackon - fixon); 
           
           objectoff = feedbackon + dfeedback - ifi;
           
           break
           

       elseif answer == 2 % if subject chose green urn 
           
           if thisurn == 1 % if thisurn is a blue urn
               
               Screen('CopyWindow', feedback_window2,window, windrect, windrect) % show "you lose" feedback window

           else % if "thisurn" is a green urn 

               Screen('CopyWindow', feedback_window1,window, windrect, windrect) % show "you win" feedback window

           end
           
           feedbackon = Screen('Flip', window, objectoff - slack);      % feedback window on 
           
           fprintf('fixation was on for %3.4f\n', feedbackon - fixon); 
           
           objectoff = feedbackon + dfeedback - ifi;
           
           break
           
       end
       
   elseif answer == 3
       
       if draw_count == drawlen % if this is the last (10th) draw ans subject requested another draw
           
           % PUT FIXATION BACK ON 
           Screen('CopyWindow', feedback_window3,window, windrect, windrect)
           feedbackon = Screen('Flip', window, objectoff - slack);      % feedback window on 
           
           fprintf('fixation was on for %3.4f\n', feedbackon - fixon); 
           
           objectoff = feedbackon + dfeedback - ifi;
           break
           
       else % if this is not the last draw and subject wants to dra again
           
           Screen('CopyWindow', fixationdisplay,window, windrect, windrect)
           objecton   = Screen('Flip', window, objectoff - slack); 
           
           fprintf('fixation was on for %3.4f\n', objecton - fixon); 
           
           objectoff = objecton + isi - ifi;
       end
   end
           
 
   % add trial related info here 
   

end % end of draw for loop

WaitSecs(2); %  1= force wait for actual pulse; 0=return this many ms after pulse
% end of block


end