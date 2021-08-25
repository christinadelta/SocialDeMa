function [set, logs] = RunTrials(set, scrn, logs, keys)

% THIS IS A SUBFUNCTION, PART OF THE "OPTIMAL STOPPING EXPERIMENTS". 

% it takes various stored information from other subfunctions to run
% one sequence

%% ---- Prepare the "global" information needed for all the tasks ---- %%

% UNPACK GLOBAL PARAMS FROM THE SETTINGS AND SCREEN STRUCTS
taskNb          = set.taskNb;       % number of task (needed to run the correct task)
total_trials    = set.trials;       % total number of trials
thistrial       = set.thistrial;    % this is the current sequence/trial number 
thisblock       = set.thisblock;     % this is the current block number 

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

%% ----- Run the correct task ------ %%

abort       = 0;
HideCursor;

% UNPACK BEADS RELATED PARAMETERS 
green       = scrn.green;
blue        = scrn.blue;
red         = scrn.red;

sequence    = set.sequence;     % this is the current sequence/trial 
thisurn     = set.urn;          % this is the current urn, contains the high prob colour beads 1=blue, 0=green
drawlen     = set.draws;        % length of each sequence (up to 10 draws)
penalty     = set.penalty;      % £0.25 loss for every draw
win         = set.win;          % £10 for every win!
balance     = set.balance;      % balance starts as zero and updates from correct/incorrect (loss) responses & draw choices

bead_dur    = set.bead_dur;     % duration of each bead in sec
response    = set.response;     % response duration in sec
dfeedback   = set.feed_dur;     % duration of feedback window in sec
confrating  = set.confrating;   % duration of the confidence rating 
fix_dur     = set.fix_dur;      % fixation duration

bluekey     = keys.code1;       % keycode for blue option
greenkey    = keys.code2;       % keycode for green option
drawkey     = keys.code3;       % keycode for draw-again option
esckey      = keys.code9;       % keycode for aborting the experiment 
    
notconfkey  = keys.code4;       % not confindent (in confidence rating)
modconfkey  = keys.code5;       % moderately confindent (in confidence rating)   
confkey     = keys.code6;       % very confindent (in confidence rating)

% ADD A FEW MORE VARIABLES
trials      = [];   % store trial info
draws       = []; % store info specifically about the sequence/draws [1-10]
draw_count  = 0;    % drawing counter
accuracy    = nan;


% SET ALL RELEVANT WINDOWS (BEAD WINDOWS, RESPONSE WINDOWS, CONFIDENCE RATING & FEEDBACK WINDOWS)

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

% CONFIDENCE RATING WINDOW
confidence_window = Screen('OpenOffscreenWindow',window);
Screen('TextSize', confidence_window, textsize);
Screen('FillRect', confidence_window, grey ,windrect);
DrawFormattedText(confidence_window, 'On a scale of 1 to 3, how confident are you for your choice?', 'center', ycenter-100, white);
DrawFormattedText(confidence_window, 'Left Arrow: Not Confident', 'center', ycenter, white);
DrawFormattedText(confidence_window, 'Down Arrow: Moderately Confident', 'center', ycenter+50, white);
DrawFormattedText(confidence_window, 'Right Arrow: Very Confident', 'center', ycenter+100, white);

% START THE TRIAL WITH A FIXATION CROSS
Screen('CopyWindow', fixationdisplay,window, windrect, windrect)
fliptime    = Screen('Flip', window); % flip fixation window
trialstart  = fliptime;

% object offset
objectoff   = trialstart + isi + randperm(jitter*1000,1)/1000 - ifi;

% 1. BEGIN DRAWING 
for thisdraw = 1:drawlen
    
    % allow subject to abort the experiment 
    [keyisdown,~,keycode] = KbCheck;
    if keyisdown && keycode(esckey)
        abort = 1;
        break;
    end
    
    % STATEMENT 1. Is it a high or a low probability draw?
    if sequence(thisdraw, 1) == 1  % high prob draw
        
        % STATEMENT 2. Is it a blue or a green bead?
        if thisurn == 1         % blue urn and blue bead
            % show blue window
            Screen('CopyWindow', blue_window, window, windrect, windrect)
            
        else % green urn and green beed
            % show green window
            Screen('CopyWindow', green_window, window, windrect, windrect)
        end
            
    elseif sequence(thisdraw, 1) == 2 % low prob draw
        
        if thisurn == 1 % blue urn and green bead
            % show green window
            Screen('CopyWindow', green_window, window, windrect, windrect)
            
        else % green urn and blue bead
            % show blue window
            Screen('CopyWindow', blue_window, window, windrect, windrect)
            
        end 
    end % end of statment 1 if 
    
    beadon      = Screen('Flip', window, objectoff - slack);                % here the current bead is fliped
    objectstart = beadon;
    
    beadstart   = objectstart - trialstart;                                 % timestamp -- start of the first bead
    beadoff     = beadon + bead_dur - ifi;                                  % bead on for 1000 ms
    
    % 2. BRING FIXATION BACK ON
    Screen('CopyWindow', fixationdisplay,window, windrect, windrect)
    fixon        = Screen('Flip', window, beadoff - slack);                 % fixation on, prepare for response prompt   
    
    fprintf('bead was on for %3.4f\n', fixon - objectstart);                % time interval from the flip of the bead until the fixation
   
    objectoff    = fixon + fix_dur - ifi;     % add jitter here?
    
    % 3. SHOW RESPONSE PROMPT 
    if thisdraw < drawlen % show response window 1
        Screen('CopyWindow', response_window1, window, windrect, windrect)
        
    else % show response window 2
        Screen('CopyWindow', response_window2, window, windrect, windrect)
      
    end % end of resp prompt if statement
    
    prompton    = Screen('Flip', window, objectoff - slack); 
    
    fprintf('fixation was on for %3.4f\n', prompton - fixon);               % time interval from the first flip until the next
    
    
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
            
            if keycode(1,bluekey) % if subject chose the blue urn 
                resp_input  = bluekey;
                rt          = secs - prompton;
                answer      = 1; % blue urn 
                respmade    = secs;
                
            elseif keycode(1,greenkey) % if subject chose the green urn 
                resp_input  = greenkey;
                rt          = secs - prompton;
                answer      = 2; % green urn
                respmade    = secs;
                
            elseif keycode(1,drawkey) % if subject chose to draw again
                resp_input  = drawkey;
                rt          = secs - prompton;
                answer      = 3; % draw-again 
                respmade    = secs;
                draw_count  = draw_count + 1;
                
            end %  
        end % if responded statement
    end % end of response while loop 
    
    promptoff   = respmade + isi - ifi;                                     % response prompt self paced or on for 2500 ms
    
   % 4. BRING BACK FIXATION. NEXT STEP:
   % A) if answer = 3, draw again,
   % B) if answer = 1 or answer = 2, ask subject to rate the confidence of their choice,
   % check if response was correct or incorrect, give appropreate feedback, break 
   % C) if draw_count == 10 and subject chose to draw again collect it as
   % wrong answer 
   
   Screen('CopyWindow', fixationdisplay,window, windrect, windrect)
   fixon        = Screen('Flip', window, promptoff - slack);                % fixation on, prepare for next draw or for feedback      
   
   fprintf('prompt was on for %3.4f\n', fixon - prompton);                  % time interval from the the flip of the prompt until the next
   
   objectoff    = fixon + fix_dur + randperm(jitter*1000,1)/1000 - ifi;     % add jitter here?
   
   % fprintf('fixation was on for %3.4f\n', objectoff); 
   
   if answer ~= 3
       
       % if subject chose an urn, before moving to feedback, ask them to
       % rate the confidence of their choice 
       Screen('CopyWindow', confidence_window, window, windrect, windrect)
       ratingon = Screen('Flip', window, objectoff - slack);  % rating window is on
       
       rat3_rt      = NaN;
       rate         = NaN;
       rate_input   = 0;
       % get response
       while rate_input == 0 && (GetSecs - ratingon) < confrating
           
           [keyisdown, secs, keycode] = KbCheck;
           pressedKeys = find(keycode);
           
           % check the response key
            if isempty(pressedKeys) % if something weird happens or subject doesn't press any of the correct keys
                rate_input  = 0; 
                rate_rt     = nan;
                rate        = nan;
                respmade    = GetSecs;

            elseif ~isempty(pressedKeys) % if subject pressed a valid key

                if keycode(1,notconfkey) % not confident
                    rate_input  = notconfkey;
                    rate_rt     = secs - ratingon;
                    rate        = 1; % not confident 
                    respmade    = secs;

                elseif keycode(1,modconfkey) % subject pressed the inanimate key
                    rate_input  = modconfkey;
                    rate_rt     = secs - ratingon;
                    rate        = 2; % moderately confident
                    respmade    = secs;

                elseif keycode(1,confkey) % subject pressed the inanimate key
                    rate_input  = confkey;
                    rate_rt     = secs - ratingon;
                    rate        = 3; % very confident
                    respmade    = secs;
                end
            end % if responded statement
        end % end of response while loop
         
       objectoff = respmade + isi + randperm(jitter*1000,1)/1000 - ifi;
       
       if answer == 1 % if subject chose blue urn
           if thisurn == 1                                                      % if thisurn is in fact a blue urn
               
               accuracy         = 1;                                                    % update the accuracy var
               
               Screen('CopyWindow', feedback_window1,window, windrect, windrect)% show "you win" feedback window

           else % if "thisurn" is a green urn 
               accuracy         = 0;
               
               Screen('CopyWindow', feedback_window2,window, windrect, windrect)% show "you lose" feedback window

           end
           
           feedbackon = Screen('Flip', window, objectoff - slack);          % feedback window on 
           
           fprintf('rating was on for %3.4f\n', feedbackon - ratingon); 
           
           objectoff = feedbackon + dfeedback - ifi;
           
           break;
           

       elseif answer == 2 % if subject chose green urn 
           
           if thisurn == 1 % if thisurn is a blue urn
               
               accuracy         = 0;
               
               Screen('CopyWindow', feedback_window2,window, windrect, windrect) % show "you lose" feedback window

           else % if "thisurn" is a green urn 

               Screen('CopyWindow', feedback_window1,window, windrect, windrect) % show "you win" feedback window
               accuracy = 1;                                                    % update the accuracy var
               
           end
           
           feedbackon = Screen('Flip', window, objectoff - slack);      % feedback window on 
           
           fprintf('rating was on for %3.4f\n', feedbackon - ratingon); 
           
           objectoff = feedbackon + dfeedback - ifi;
           
           break
           
       end
       
   elseif answer == 3
       
       if draw_count == drawlen % if this is the last (10th) draw ans subject requested another draw
           
           accuracy     = 0;
           % PUT FIXATION BACK ON 
           Screen('CopyWindow', feedback_window3,window, windrect, windrect)
           feedbackon   = Screen('Flip', window, objectoff - slack);      % feedback window on 
           
           fprintf('fixation was on for %3.4f\n', feedbackon - fixon); 
           
           objectoff    = feedbackon + dfeedback - ifi;
           break
           
       else % if this is not the last draw and subject wants to draw again
           
           Screen('CopyWindow', fixationdisplay,window, windrect, windrect)
           objecton     = Screen('Flip', window, objectoff - slack); 
           
           fprintf('fixation was on for %3.4f\n', objecton - fixon); 
           
           objectoff    = objecton + isi - ifi;
           
       end
   end
           
   % add sequence/draws related info here 
   draws(thisdraw).session      = thisession;
   draws(thisdraw).block        = thisblock;
   draws(thisdraw).trialnumber  = thistrial;
   draws(thisdraw).trialonset   = trialstart;
   draws(thisdraw).currentdraw  = thisdraw;
   draws(thisdraw).rt           = rt;
   
   if abort; fclose('all');break; end 

end % end of draw for loop

% UPDATE BALANCE 
% 1. First subtract the penalty (0.25p for every draw)
balance = balance - (penalty * draw_count);

% 2. Now add winnings of sabtract loss based on the accuracy and loss function 
if accuracy == 1
    balance = balance + win;
    
elseif accuracy == 0 && unique(sequence(:,2)) == 0
    balance = balance;
    
elseif accuracy == 0 && unique(sequence(:,2)) == 1
    balance = balance - win;
    
end

% update the current balance
set.balance                  = balance; 
   

% add trial related info here 
trials.session      = thisession;
trials.block        = thisblock;
trials.trialnumber  = thistrial;
trials.trialonset   = trialstart;
trials.urntype      = thisurn;
trials.sequence     = sequence(:,1)';
trials.loss         = unique(sequence(:,2));
trials.draws        = draw_count;
trials.response     = answer;
trials.accuracy     = accuracy;
trials.balance      = balance;

WaitSecs(2); %  1= force wait for actual pulse; 0=return this many ms after pulse

% store draws and trials info in the logs mat file 
logs.trialstart         = trialstart;
logs.draws              = draws;
set.trials              = trials;

sub_drawslog             = fullfile(resfolder,sprintf(logs.drawslog,sub,taskname,thisblock,thistrial,thisession));
save(sub_drawslog,'logs');


% end of block


end