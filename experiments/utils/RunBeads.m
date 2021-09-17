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
grey            = scrn.grey;
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

bluekey     = set.code1;       % keycode for blue option
greenkey    = set.code2;       % keycode for green option
drawkey     = set.code3;       % keycode for draw-again option
esckey      = set.code21;      % keycode for aborting the experiment 
    
notconfkey  = set.code4;       % not confindent (in confidence rating)
modconfkey  = set.code5;       % moderately confindent (in confidence rating)   
confkey     = set.code6;       % very confindent (in confidence rating)

% ADD A FEW MORE VARIABLES
trials      = []; % store trial info
draws       = []; % store info specifically about the sequence/draws [1-10]
previous    = []; % store the previous draws here to show at the bottom of the screen
colours     = []; % store the previous draw-colours here to show at the bottom of the screen
draw_count  = 0;  % drawing counter
accuracy    = nan;

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
    
    xcenters        = xcenter - 600; % start adding the darws prices at xcenter - 600
    
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
            DrawFormattedText(window, 'Previous Draws', xcenter - 600, ycenter+250, white); 
            DrawFormattedText(window, previous{i}, xcenters, ycenter+350, colours{i}); 

            % update xcenters, so that previous prices are not
            % displayed on top of each other
            xcenters = xcenters+80;
            
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
    fixation_onset        = Screen('Flip', window, bead_offset - slack);        % fixation on, prepare for response prompt   
    
    object_offset    = fixation_onset + fix_dur - ifi;     % add jitter here?
    
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
           
           Screen('CopyWindow', fixationdisplay,window, windrect, windrect)
           object_onset     = Screen('Flip', window, object_offset - slack); 

           object_offset    = object_onset + isi + randperm(jitter*1000,1)/1000 - ifi;% add jitter here?
           
       else
       
           % if subject chose an urn, before moving to feedback, ask them to
           % rate the confidence of their choice 
           Screen('CopyWindow', confidence_window, window, windrect, windrect)
           ratingon = Screen('Flip', window, object_offset - slack);  % rating window is on

           % send confidence screen trigger
            if EEG == 1 
                sp.sendTrigger(trigger9)
            end

           rat3_rt              = NaN;
           rate                 = NaN;
           rate_input           = 0;
           responseTrigNotSent  = 1;

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

                        % send not confident trigger
                        if EEG == 1 && responseTrigNotSent==1
                            sp.sendTrigger(trigger10);
                            responseTrigNotSent=0;
                        end

                    elseif keycode(1,modconfkey) % subject pressed the inanimate key
                        rate_input  = modconfkey;
                        rate_rt     = secs - ratingon;
                        rate        = 2; % moderately confident
                        respmade    = secs;

                        % send moderately confident trigger
                        if EEG == 1 && responseTrigNotSent==1
                            sp.sendTrigger(trigger11);
                            responseTrigNotSent=0;
                        end

                    elseif keycode(1,confkey) % subject pressed the inanimate key
                        rate_input  = confkey;
                        rate_rt     = secs - ratingon;
                        rate        = 3; % very confident
                        respmade    = secs;

                        % send very confident trigger
                        if EEG == 1 && responseTrigNotSent==1
                            sp.sendTrigger(trigger12);
                            responseTrigNotSent=0;
                        end

                    end
                end % if responded statement
            end % end of response while loop
            
            object_offset = respmade + isi + randperm(jitter*1000,1)/1000 - ifi; %self-paced
           
            % ADD A FIXATION BEFORE FEEDBACK
            Screen('CopyWindow', fixationdisplay,window, windrect, windrect)
            object_onset     = Screen('Flip', window, object_offset - slack); 

            object_offset    = object_onset + fix_dur - ifi;
           

           if answer == 1 % if subject chose blue urn
               if thisurn == 1                                                      % if thisurn is in fact a blue urn

                   accuracy         = 1;                                            % update the accuracy var

                   Screen('CopyWindow', feedback_window1,window, windrect, windrect)% show "you win" feedback window
                   feedback_onset = Screen('Flip', window, object_offset - slack);  % feedback window on 

                   % send you win trigger
                    if EEG == 1
                        sp.sendTrigger(trigger13);
                    end

               else % if "thisurn" is a green urn 
                   accuracy         = 0;

                   Screen('CopyWindow', feedback_window2,window, windrect, windrect) % show "you lose" feedback window
                   feedback_onset = Screen('Flip', window, object_offset - slack);   % feedback window on 

                   % send you lose trigger
                    if EEG == 1
                        sp.sendTrigger(trigger14);
                    end
               end
               object_offset = feedback_onset + dfeedback - ifi;

               break;


           elseif answer == 2 % if subject chose green urn 

               if thisurn == 1 % if thisurn is a blue urn

                   accuracy         = 0;

                   Screen('CopyWindow', feedback_window2,window, windrect, windrect) % show "you lose" feedback window
                   feedback_onset = Screen('Flip', window, object_offset - slack);   % feedback window on 
                   % send you lose trigger
                    if EEG == 1
                        sp.sendTrigger(trigger14);
                    end

               else % if "thisurn" is a green urn 

                   Screen('CopyWindow', feedback_window1,window, windrect, windrect) % show "you win" feedback window
                   accuracy = 1;                                                     % update the accuracy var
                   feedback_onset = Screen('Flip', window, object_offset - slack);   % feedback window on 
                   % send you lose trigger
                   if EEG == 1
                       sp.sendTrigger(trigger13);
                   end
               end
               object_offset = feedback_onset + dfeedback - ifi;
               break

           end
       end
       
   elseif answer == 3
       
       if draw_count == drawlen % if this is the last (10th) draw ans subject requested another draw
           
           accuracy     = 0;
           % PUT FIXATION BACK ON 
           Screen('CopyWindow', feedback_window3,window, windrect, windrect)
           feedback_onset   = Screen('Flip', window, object_offset - slack); % feedback window on 
           % send you lose (out of draws) trigger
           if EEG == 1
               sp.sendTrigger(trigger15);
           end
           object_offset    = feedback_onset + dfeedback - ifi;
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

% 2. Now add winnings of sabtract loss based on the accuracy and loss function 
if accuracy == 1
    balance = balance + win;
  
elseif accuracy == 0 
    balance = balance - win;
    
end

% update the current balance
set.balance         = balance; 
   
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